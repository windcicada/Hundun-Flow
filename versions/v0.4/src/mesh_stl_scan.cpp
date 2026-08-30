// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_mesh.hpp"

#include "mesh_stl_scan_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <ios>
#include <limits>
#include <locale>
#include <new>
#include <cstdlib>
#include <string>
#include <string_view>
#include <sys/stat.h>
#include <thread>
#include <type_traits>
#include <unistd.h>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
std::atomic<std::uint64_t> g_stl_open_count{0U};
std::atomic<std::size_t> g_thread_launches_before_failure{
    std::numeric_limits<std::size_t>::max()};
#endif

constexpr std::uint64_t kFnvOffset = UINT64_C(14695981039346656037);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);
constexpr std::uint64_t kMaxBinReferences = UINT64_C(32) * 1024U * 1024U;
constexpr double kRoundoffScale = 64.0;

Status invalid_plan(detail::StlScanDetail value) noexcept {
  return {StatusCode::invalid_plan, static_cast<std::uint32_t>(value)};
}

Status io_failure(detail::StlScanDetail value) noexcept {
  return {StatusCode::io_failure, static_cast<std::uint32_t>(value)};
}

Status allocation_failure() noexcept {
  return {StatusCode::allocation_failure,
          static_cast<std::uint32_t>(detail::stl_detail_allocation)};
}

Status mpi_failure() noexcept {
  return {StatusCode::mpi_failure,
          static_cast<std::uint32_t>(detail::stl_detail_collective)};
}

class Hash64 {
 public:
  void byte(std::uint8_t value) noexcept {
    value_ ^= value;
    value_ *= kFnvPrime;
  }

  template <class Integer>
  void integer(Integer value) noexcept {
    using Unsigned = std::make_unsigned_t<Integer>;
    const Unsigned bits = static_cast<Unsigned>(value);
    for (std::size_t shift = 0U; shift < sizeof(bits) * 8U; shift += 8U) {
      byte(static_cast<std::uint8_t>((bits >> shift) & Unsigned{0xffU}));
    }
  }

  void real(double value) noexcept {
    std::uint64_t bits{};
    static_assert(sizeof(bits) == sizeof(value));
    std::memcpy(&bits, &value, sizeof(bits));
    integer(bits);
  }

  PlanFingerprint finish() const noexcept {
    return value_ == 0U ? PlanFingerprint{1U} : value_;
  }

 private:
  std::uint64_t value_{kFnvOffset};
};

int axis_index(CartesianAxis axis) noexcept {
  switch (axis) {
    case CartesianAxis::x:
      return 0;
    case CartesianAxis::y:
      return 1;
    case CartesianAxis::z:
      return 2;
  }
  return -1;
}

std::array<int, 2U> transverse_axes(int axis) noexcept {
  if (axis == 0) {
    return {1, 2};
  }
  if (axis == 1) {
    return {0, 2};
  }
  return {0, 1};
}

std::int32_t component(Int3 value, int axis) noexcept {
  return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

double component(Real3 value, int axis) noexcept {
  return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

bool checked_product(std::size_t a, std::size_t b,
                     std::size_t& out) noexcept {
  if (a != 0U && b > std::numeric_limits<std::size_t>::max() / a) {
    return false;
  }
  out = a * b;
  return true;
}

bool checked_u64_add(std::uint64_t& value, std::uint64_t amount) noexcept {
  if (amount > std::numeric_limits<std::uint64_t>::max() - value) {
    return false;
  }
  value += amount;
  return true;
}

bool checked_u64_product(std::uint64_t a, std::uint64_t b,
                         std::uint64_t& out) noexcept {
  if (a != 0U && b > std::numeric_limits<std::uint64_t>::max() / a) {
    return false;
  }
  out = a * b;
  return true;
}

bool valid_budget(StlScanBudget budget) noexcept {
  return budget.max_persistent_bytes_per_rank != 0U &&
         budget.max_peak_bytes_per_rank != 0U &&
         budget.max_peak_bytes_per_rank >=
             budget.max_persistent_bytes_per_rank &&
         budget.max_bin_references != 0U &&
         budget.max_events_per_line != 0U &&
         budget.worker_threads != 0U && budget.worker_threads <= 1024U &&
         budget.max_bin_references <=
             static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()) &&
         budget.max_events_per_line <=
             static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max());
}

bool inject_thread_launch_failure() noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  std::size_t remaining =
      g_thread_launches_before_failure.load(std::memory_order_relaxed);
  for (;;) {
    if (remaining == std::numeric_limits<std::size_t>::max()) {
      return false;
    }
    if (remaining == 0U) {
      return true;
    }
    if (g_thread_launches_before_failure.compare_exchange_weak(
            remaining, remaining - 1U, std::memory_order_relaxed,
            std::memory_order_relaxed)) {
      return false;
    }
  }
#else
  return false;
#endif
}

std::uint64_t bin_reference_limit(StlScanBudget budget) noexcept {
  return std::min(kMaxBinReferences, budget.max_bin_references);
}

Status initial_memory_gate(std::size_t triangles, std::size_t lines,
                           std::uint64_t source_bytes,
                           StlScanBudget budget) noexcept {
  if (!valid_budget(budget)) {
    return invalid_plan(detail::stl_detail_budget);
  }
  const std::uint64_t triangle_count = static_cast<std::uint64_t>(triangles);
  const std::uint64_t line_count = static_cast<std::uint64_t>(lines);
  std::uint64_t soa_bytes{};
  std::uint64_t line_bytes{};
  std::uint64_t triangle_input_bytes{};
  std::uint64_t coordinate_wire_bytes{};
  std::uint64_t triangle_hash_bytes{};
  std::uint64_t event_bytes{};
  std::uint64_t worker_output_bytes{};
  std::uint64_t worker_control_bytes{};
  std::uint64_t maximum_span_bytes{};
  const std::uint64_t reference_limit = bin_reference_limit(budget);
  // TriangleSoA owns 24 double arrays.  Line offsets and at most the configured
  // total span-bound count are persistent upper bounds.
  if (!checked_u64_product(triangle_count, UINT64_C(24) * sizeof(double),
                           soa_bytes) ||
      !checked_u64_product(line_count + 1U, sizeof(std::size_t), line_bytes) ||
      !checked_u64_product(triangle_count, sizeof(TriangleInput),
                           triangle_input_bytes) ||
      !checked_u64_product(triangle_count, UINT64_C(9) * sizeof(double),
                           coordinate_wire_bytes) ||
      !checked_u64_product(triangle_count, sizeof(std::uint64_t),
                           triangle_hash_bytes) ||
      // Event is one binary64 coordinate plus sign/key and ABI padding; use a
      // conservative 32-byte planning charge independent of host layout.
      !checked_u64_product(budget.max_events_per_line, UINT64_C(32),
                           event_bytes) ||
      !checked_u64_product(event_bytes, budget.worker_threads, event_bytes) ||
      !checked_u64_product(reference_limit, sizeof(double),
                           maximum_span_bytes) ||
      !checked_u64_product(line_count, sizeof(std::size_t),
                           worker_output_bytes) ||
      !checked_u64_add(worker_output_bytes, maximum_span_bytes) ||
      !checked_u64_product(budget.worker_threads, UINT64_C(256),
                           worker_control_bytes)) {
    return invalid_plan(detail::stl_detail_budget);
  }
  std::uint64_t persistent = 0U;
  if (!checked_u64_add(persistent, soa_bytes) ||
      !checked_u64_add(persistent, line_bytes) ||
      !checked_u64_add(persistent, maximum_span_bytes) ||
      persistent > budget.max_persistent_bytes_per_rank) {
    return invalid_plan(detail::stl_detail_budget);
  }
  std::uint64_t peak = persistent;
  std::uint64_t bin_offsets_bytes{};
  std::uint64_t bin_reference_bytes{};
  if (!checked_u64_product(line_count + 1U, sizeof(std::size_t),
                           bin_offsets_bytes) ||
      !checked_u64_product(reference_limit, sizeof(std::size_t),
                           bin_reference_bytes) ||
      !checked_u64_add(peak, triangle_input_bytes) ||
      !checked_u64_add(peak, coordinate_wire_bytes) ||
      !checked_u64_add(peak, source_bytes) ||
      !checked_u64_add(peak, triangle_hash_bytes) ||
      !checked_u64_add(peak, bin_offsets_bytes) ||
      !checked_u64_add(peak, bin_offsets_bytes) ||
      !checked_u64_add(peak, bin_reference_bytes) ||
      !checked_u64_add(peak, event_bytes) ||
      !checked_u64_add(peak, worker_output_bytes) ||
      !checked_u64_add(peak, worker_control_bytes) ||
      // During deterministic merge, worker-owned spans and the final plan's
      // span storage coexist.  Charge another complete span upper bound so a
      // skewed chunk cannot cross the declared peak between reserve and move.
      !checked_u64_add(peak, maximum_span_bytes) ||
      peak > budget.max_peak_bytes_per_rank) {
    return invalid_plan(detail::stl_detail_budget);
  }
  return {};
}

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

Real3 subtract(Real3 lhs, Real3 rhs) noexcept {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

Real3 cross(Real3 lhs, Real3 rhs) noexcept {
  return {lhs.y * rhs.z - lhs.z * rhs.y,
          lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

double dot(Real3 lhs, Real3 rhs) noexcept {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

double norm_squared(Real3 value) noexcept { return dot(value, value); }

bool same_coordinate(double lhs, double rhs) noexcept {
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs - rhs) <=
         kRoundoffScale * std::numeric_limits<double>::epsilon() * scale;
}

struct Event {
  double coordinate{};
  std::int8_t normal_sign{};
  std::uint64_t triangle_key{};
};

bool event_less(const Event& lhs, const Event& rhs) noexcept {
  if (lhs.coordinate != rhs.coordinate) {
    return lhs.coordinate < rhs.coordinate;
  }
  if (lhs.normal_sign != rhs.normal_sign) {
    return lhs.normal_sign < rhs.normal_sign;
  }
  return lhs.triangle_key < rhs.triangle_key;
}

bool coordinate_bits_less(Real3 lhs, Real3 rhs) noexcept {
  if (lhs.x != rhs.x) {
    return lhs.x < rhs.x;
  }
  if (lhs.y != rhs.y) {
    return lhs.y < rhs.y;
  }
  return lhs.z < rhs.z;
}

std::array<Real3, 3U> canonical_vertices(const TriangleInput& triangle) noexcept {
  std::array<Real3, 3U> vertices{triangle.a, triangle.b, triangle.c};
  std::sort(vertices.begin(), vertices.end(), coordinate_bits_less);
  return vertices;
}

std::uint64_t triangle_key(Real3 vertex, Real3 edge1, Real3 edge2) noexcept {
  Hash64 hash;
  hash.real(vertex.x);
  hash.real(vertex.y);
  hash.real(vertex.z);
  hash.real(edge1.x);
  hash.real(edge1.y);
  hash.real(edge1.z);
  hash.real(edge2.x);
  hash.real(edge2.y);
  hash.real(edge2.z);
  return hash.finish();
}

std::uint32_t read_u32(const std::uint8_t* bytes) noexcept {
  return static_cast<std::uint32_t>(bytes[0]) |
         (static_cast<std::uint32_t>(bytes[1]) << 8U) |
         (static_cast<std::uint32_t>(bytes[2]) << 16U) |
         (static_cast<std::uint32_t>(bytes[3]) << 24U);
}

float read_f32(const std::uint8_t* bytes) noexcept {
  const std::uint32_t bits = read_u32(bytes);
  float value{};
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

class AsciiCursor {
 public:
  explicit AsciiCursor(Span<const std::uint8_t> bytes) noexcept
      : current_(reinterpret_cast<const char*>(bytes.data)),
        end_(current_ + bytes.size) {}

  bool token(std::string_view& out) noexcept {
    skip_space();
    const char* begin = current_;
    while (current_ != end_ && !space(*current_)) {
      ++current_;
    }
    out = {begin, static_cast<std::size_t>(current_ - begin)};
    return !out.empty();
  }

  void skip_line() noexcept {
    while (current_ != end_ && *current_ != '\n') {
      ++current_;
    }
    if (current_ != end_) {
      ++current_;
    }
  }

  bool finished() noexcept {
    skip_space();
    return current_ == end_;
  }

 private:
  static bool space(char value) noexcept {
    return value == ' ' || value == '\t' || value == '\r' || value == '\n';
  }
  void skip_space() noexcept {
    while (current_ != end_ && space(*current_)) {
      ++current_;
    }
  }
  const char* current_{};
  const char* end_{};
};

bool parse_real_token(std::string_view token, double& out) noexcept {
  if (token.empty() || token.size() > 128U) {
    return false;
  }
  std::array<char, 129U> buffer{};
  std::copy(token.begin(), token.end(), buffer.begin());
  char* parsed_end = nullptr;
  errno = 0;
  out = std::strtod(buffer.data(), &parsed_end);
  return errno != ERANGE && parsed_end == buffer.data() + token.size() &&
         std::isfinite(out);
}

bool read_real(AsciiCursor& cursor, double& out) noexcept {
  std::string_view token;
  return cursor.token(token) && parse_real_token(token, out);
}

bool expect_token(AsciiCursor& cursor, std::string_view expected) noexcept {
  std::string_view token;
  return cursor.token(token) && token == expected;
}

bool parse_ascii(Span<const std::uint8_t> bytes,
                 std::vector<TriangleInput>& triangles) {
  if (bytes.data == nullptr || bytes.size == 0U ||
      std::find(bytes.data, bytes.data + bytes.size, std::uint8_t{0U}) !=
          bytes.data + bytes.size) {
    return false;
  }

  AsciiCursor cursor(bytes);
  std::string_view word;
  if (!cursor.token(word) || word != "solid") {
    return false;
  }
  cursor.skip_line();

  triangles.clear();
  for (;;) {
    if (!cursor.token(word)) {
      return false;
    }
    if (word == "endsolid") {
      cursor.skip_line();
      return cursor.finished() && !triangles.empty();
    }
    if (word != "facet" || !expect_token(cursor, "normal")) {
      return false;
    }
    Real3 file_normal{};
    if (!read_real(cursor, file_normal.x) ||
        !read_real(cursor, file_normal.y) ||
        !read_real(cursor, file_normal.z) || !finite(file_normal)) {
      return false;
    }
    if (!expect_token(cursor, "outer") || !expect_token(cursor, "loop")) {
      return false;
    }
    TriangleInput triangle{};
    Real3* vertices[3U]{&triangle.a, &triangle.b, &triangle.c};
    for (Real3* vertex : vertices) {
      if (!expect_token(cursor, "vertex") ||
          !read_real(cursor, vertex->x) || !read_real(cursor, vertex->y) ||
          !read_real(cursor, vertex->z) || !finite(*vertex)) {
        return false;
      }
    }
    if (!expect_token(cursor, "endloop") ||
        !expect_token(cursor, "endfacet")) {
      return false;
    }
    if (triangles.size() >= detail::kMaxStlTriangles) {
      return false;
    }
    triangles.push_back(triangle);
  }
}

bool appears_text_stl(Span<const std::uint8_t> bytes) noexcept {
  if (bytes.size < 5U || bytes.data == nullptr ||
      std::memcmp(bytes.data, "solid", 5U) != 0) {
    return false;
  }
  for (std::size_t index = 0U; index < bytes.size; ++index) {
    const std::uint8_t value = bytes.data[index];
    if (value == 0U || (value < 0x09U) ||
        (value > 0x0dU && value < 0x20U)) {
      return false;
    }
  }
  return true;
}

bool parse_binary(Span<const std::uint8_t> bytes, std::uint32_t count,
                  std::vector<TriangleInput>& triangles) {
  constexpr std::uint64_t kHeaderBytes = 84U;
  constexpr std::uint64_t kRecordBytes = 50U;
  const std::uint64_t expected =
      kHeaderBytes + static_cast<std::uint64_t>(count) * kRecordBytes;
  if (count == 0U || count > detail::kMaxStlTriangles ||
      expected != bytes.size) {
    return false;
  }
  triangles.clear();
  triangles.reserve(count);
  for (std::size_t triangle_index = 0U; triangle_index < count;
       ++triangle_index) {
    const std::uint8_t* record =
        bytes.data + 84U + triangle_index * 50U;
    const Real3 file_normal{static_cast<double>(read_f32(record)),
                            static_cast<double>(read_f32(record + 4U)),
                            static_cast<double>(read_f32(record + 8U))};
    if (!finite(file_normal)) {
      return false;
    }
    TriangleInput triangle{};
    Real3* vertices[3U]{&triangle.a, &triangle.b, &triangle.c};
    for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
      const std::uint8_t* source = record + 12U + vertex * 12U;
      *vertices[vertex] = {static_cast<double>(read_f32(source)),
                           static_cast<double>(read_f32(source + 4U)),
                           static_cast<double>(read_f32(source + 8U))};
      if (!finite(*vertices[vertex])) {
        return false;
      }
    }
    triangles.push_back(triangle);
  }
  return true;
}

Status read_stl_file(const std::filesystem::path& case_root,
                     const std::filesystem::path& relative,
                     std::size_t line_count,
                     StlScanBudget budget,
                     std::vector<std::uint8_t>& bytes,
                     std::size_t& triangle_upper_bound) {
  if (relative.empty() || relative.is_absolute() || relative.has_parent_path() ||
      relative.filename() != relative || relative == "." || relative == "..") {
    return io_failure(detail::stl_detail_path);
  }

  const int root = ::open(case_root.c_str(),
                          O_RDONLY | O_CLOEXEC | O_DIRECTORY | O_NOFOLLOW);
  if (root < 0) {
    return io_failure(detail::stl_detail_open);
  }
  const int descriptor = ::openat(
      root, relative.c_str(),
      O_RDONLY | O_CLOEXEC | O_NOFOLLOW | O_NONBLOCK | O_NOCTTY);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  g_stl_open_count.fetch_add(1U, std::memory_order_relaxed);
#endif
  const int open_error = errno;
  static_cast<void>(::close(root));
  if (descriptor < 0) {
    static_cast<void>(open_error);
    return io_failure(detail::stl_detail_open);
  }

  struct stat metadata {};
  if (::fstat(descriptor, &metadata) != 0 || !S_ISREG(metadata.st_mode) ||
      metadata.st_size <= 0 ||
      static_cast<std::uint64_t>(metadata.st_size) > detail::kMaxStlBytes ||
      static_cast<std::uint64_t>(metadata.st_size) >
          budget.max_peak_bytes_per_rank) {
    static_cast<void>(::close(descriptor));
    return io_failure(detail::stl_detail_size);
  }
  const std::uint64_t source_bytes =
      static_cast<std::uint64_t>(metadata.st_size);
  // A binary facet is exactly 50 bytes after its header; a valid ASCII facet
  // needs more bytes.  ceil(file_size/50) is therefore a conservative format-
  // independent bound available before reading or parsing the payload.
  const std::uint64_t triangle_bound = source_bytes / UINT64_C(50) + 1U;
  if (triangle_bound > detail::kMaxStlTriangles ||
      triangle_bound >
          static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    static_cast<void>(::close(descriptor));
    return invalid_plan(detail::stl_detail_budget);
  }
  const Status preflight = initial_memory_gate(
      static_cast<std::size_t>(triangle_bound), line_count, source_bytes,
      budget);
  if (!preflight) {
    static_cast<void>(::close(descriptor));
    return preflight;
  }
  triangle_upper_bound = static_cast<std::size_t>(triangle_bound);
  bytes.resize(static_cast<std::size_t>(metadata.st_size));
  std::size_t total = 0U;
  while (total < bytes.size()) {
    const ssize_t count =
        ::read(descriptor, bytes.data() + total, bytes.size() - total);
    if (count < 0 && errno == EINTR) {
      continue;
    }
    if (count <= 0) {
      static_cast<void>(::close(descriptor));
      return io_failure(detail::stl_detail_read);
    }
    total += static_cast<std::size_t>(count);
  }
  std::uint8_t extra{};
  ssize_t extra_count{};
  do {
    extra_count = ::read(descriptor, &extra, 1U);
  } while (extra_count < 0 && errno == EINTR);
  static_cast<void>(::close(descriptor));
  return extra_count == 0 ? Status{} : io_failure(detail::stl_detail_read);
}

Status broadcast_doubles(MPI_Comm communicator, int root, double* data,
                         std::size_t count) noexcept {
  double* position = data;
  std::size_t remaining = count;
  while (remaining != 0U) {
    const int chunk = static_cast<int>(std::min<std::size_t>(
        remaining, static_cast<std::size_t>(std::numeric_limits<int>::max())));
    if (MPI_Bcast(position, chunk, MPI_DOUBLE, root, communicator) !=
        MPI_SUCCESS) {
      return mpi_failure();
    }
    position += chunk;
    remaining -= static_cast<std::size_t>(chunk);
  }
  return {};
}

Status collective_status(MPI_Comm communicator, int rank,
                         Status local) noexcept {
  int selected = local ? std::numeric_limits<int>::max() : rank;
  int first_failure = selected;
  if (MPI_Allreduce(&selected, &first_failure, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return mpi_failure();
  }
  if (first_failure == std::numeric_limits<int>::max()) {
    return {};
  }
  std::array<std::uint32_t, 2U> wire{};
  if (rank == first_failure) {
    wire = {static_cast<std::uint32_t>(local.code), local.detail};
  }
  if (MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_UINT32_T,
                first_failure, communicator) != MPI_SUCCESS) {
    return mpi_failure();
  }
  return {static_cast<StatusCode>(wire[0]), wire[1]};
}

}  // namespace

namespace detail {

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
void reset_stl_open_count_for_test() noexcept {
  g_stl_open_count.store(0U, std::memory_order_relaxed);
}

std::uint64_t stl_open_count_for_test() noexcept {
  return g_stl_open_count.load(std::memory_order_relaxed);
}

void fail_thread_launch_after_for_test(
    std::size_t successful_launches) noexcept {
  g_thread_launches_before_failure.store(successful_launches,
                                         std::memory_order_relaxed);
}

void reset_thread_launch_failure_for_test() noexcept {
  g_thread_launches_before_failure.store(
      std::numeric_limits<std::size_t>::max(), std::memory_order_relaxed);
}
#endif

LineChunk fixed_line_chunk(std::size_t line_count, std::size_t worker,
                           std::size_t worker_count) noexcept {
  if (worker_count == 0U || worker >= worker_count) {
    return {};
  }
  const std::size_t quotient = line_count / worker_count;
  const std::size_t remainder = line_count % worker_count;
  const std::size_t begin =
      worker * quotient + std::min(worker, remainder);
  return {begin, begin + quotient + (worker < remainder ? 1U : 0U)};
}

Status parse_stl_bytes(Span<const std::uint8_t> bytes,
                       std::vector<TriangleInput>& triangles) noexcept {
  try {
    if (bytes.data == nullptr || bytes.size == 0U ||
        bytes.size > kMaxStlBytes) {
      return io_failure(stl_detail_size);
    }
    bool parsed = false;
    if (bytes.size >= 84U) {
      const std::uint32_t count = read_u32(bytes.data + 80U);
      const std::uint64_t expected =
          UINT64_C(84) + static_cast<std::uint64_t>(count) * UINT64_C(50);
      if (expected == bytes.size || !appears_text_stl(bytes)) {
        parsed = parse_binary(bytes, count, triangles);
      } else {
        parsed = parse_ascii(bytes, triangles);
      }
    } else {
      parsed = parse_ascii(bytes, triangles);
    }
    if (!parsed) {
      triangles.clear();
      return io_failure(stl_detail_syntax);
    }
    return {};
  } catch (const std::bad_alloc&) {
    triangles.clear();
    return allocation_failure();
  } catch (...) {
    triangles.clear();
    return io_failure(stl_detail_syntax);
  }
}

}  // namespace detail

Status StlScanCompiler::compile_triangles(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch,
    Span<const TriangleInput> triangles, CartesianAxis scan_axis,
    StlScanBudget budget,
    StlScanPlan& out) noexcept {
  try {
    const int scan = axis_index(scan_axis);
    const Int3 global = geometry.global_cells();
    if (scan < 0 || geometry.fingerprint() == 0U || global.x <= 0 ||
        global.y <= 0 || global.z <= 0) {
      return invalid_plan(detail::stl_detail_geometry);
    }
    if (triangles.size != 0U && triangles.data == nullptr) {
      return invalid_plan(detail::stl_detail_triangle);
    }
    for (int axis = 0; axis < 3; ++axis) {
      const std::int32_t begin = component(patch.begin, axis);
      const std::int32_t count = component(patch.cells, axis);
      const std::int32_t extent = component(global, axis);
      if (begin < 0 || count <= 0 || begin > extent || count > extent - begin ||
          geometry.axis(static_cast<CartesianAxis>(axis)).centres().size !=
              static_cast<std::size_t>(extent)) {
        return invalid_plan(detail::stl_detail_patch);
      }
    }

    const auto transverse = transverse_axes(scan);
    const std::size_t first_count = static_cast<std::size_t>(
        component(patch.cells, transverse[0]));
    const std::size_t second_count = static_cast<std::size_t>(
        component(patch.cells, transverse[1]));
    std::size_t line_count{};
    if (!checked_product(first_count, second_count, line_count)) {
      return invalid_plan(detail::stl_detail_patch);
    }
    const std::size_t triangle_count = triangles.size;
    if (triangle_count > detail::kMaxStlTriangles) {
      return invalid_plan(detail::stl_detail_triangle);
    }
    const Status memory_gate =
        initial_memory_gate(triangle_count, line_count, 0U, budget);
    if (!memory_gate) {
      return memory_gate;
    }

    StlScanPlan candidate;
    candidate.scan_axis_ = scan_axis;
    candidate.patch_begin_ = patch.begin;
    candidate.patch_cells_ = patch.cells;
    candidate.geometry_fingerprint_ = geometry.fingerprint();

    TriangleSoA& soa = candidate.triangles_;
    soa.ax_.reserve(triangle_count);
    soa.ay_.reserve(triangle_count);
    soa.az_.reserve(triangle_count);
    soa.bx_.reserve(triangle_count);
    soa.by_.reserve(triangle_count);
    soa.bz_.reserve(triangle_count);
    soa.cx_.reserve(triangle_count);
    soa.cy_.reserve(triangle_count);
    soa.cz_.reserve(triangle_count);
    soa.e1x_.reserve(triangle_count);
    soa.e1y_.reserve(triangle_count);
    soa.e1z_.reserve(triangle_count);
    soa.e2x_.reserve(triangle_count);
    soa.e2y_.reserve(triangle_count);
    soa.e2z_.reserve(triangle_count);
    soa.nx_.reserve(triangle_count);
    soa.ny_.reserve(triangle_count);
    soa.nz_.reserve(triangle_count);
    soa.min_x_.reserve(triangle_count);
    soa.min_y_.reserve(triangle_count);
    soa.min_z_.reserve(triangle_count);
    soa.max_x_.reserve(triangle_count);
    soa.max_y_.reserve(triangle_count);
    soa.max_z_.reserve(triangle_count);

    std::vector<std::uint64_t> canonical_triangle_hashes;
    canonical_triangle_hashes.reserve(triangle_count);
    for (std::size_t index = 0U; index < triangle_count; ++index) {
      const TriangleInput triangle = triangles.data[index];
      if (!finite(triangle.a) || !finite(triangle.b) || !finite(triangle.c)) {
        return invalid_plan(detail::stl_detail_triangle);
      }
      const Real3 edge1 = subtract(triangle.b, triangle.a);
      const Real3 edge2 = subtract(triangle.c, triangle.a);
      const Real3 edge3 = subtract(triangle.c, triangle.b);
      const Real3 raw_normal = cross(edge1, edge2);
      const double edge_scale_squared =
          std::max({norm_squared(edge1), norm_squared(edge2),
                    norm_squared(edge3)});
      const double normal_magnitude = std::sqrt(norm_squared(raw_normal));
      const double minimum_normal =
          kRoundoffScale * std::numeric_limits<double>::epsilon() *
          edge_scale_squared;
      if (!(edge_scale_squared > 0.0) || !std::isfinite(edge_scale_squared) ||
          !std::isfinite(normal_magnitude) ||
          !(normal_magnitude > minimum_normal)) {
        return invalid_plan(detail::stl_detail_triangle);
      }
      const Real3 normal{raw_normal.x / normal_magnitude,
                         raw_normal.y / normal_magnitude,
                         raw_normal.z / normal_magnitude};
      const Real3 minimum{std::min({triangle.a.x, triangle.b.x, triangle.c.x}),
                          std::min({triangle.a.y, triangle.b.y, triangle.c.y}),
                          std::min({triangle.a.z, triangle.b.z, triangle.c.z})};
      const Real3 maximum{std::max({triangle.a.x, triangle.b.x, triangle.c.x}),
                          std::max({triangle.a.y, triangle.b.y, triangle.c.y}),
                          std::max({triangle.a.z, triangle.b.z, triangle.c.z})};

      soa.ax_.push_back(triangle.a.x);
      soa.ay_.push_back(triangle.a.y);
      soa.az_.push_back(triangle.a.z);
      soa.bx_.push_back(triangle.b.x);
      soa.by_.push_back(triangle.b.y);
      soa.bz_.push_back(triangle.b.z);
      soa.cx_.push_back(triangle.c.x);
      soa.cy_.push_back(triangle.c.y);
      soa.cz_.push_back(triangle.c.z);
      soa.e1x_.push_back(edge1.x);
      soa.e1y_.push_back(edge1.y);
      soa.e1z_.push_back(edge1.z);
      soa.e2x_.push_back(edge2.x);
      soa.e2y_.push_back(edge2.y);
      soa.e2z_.push_back(edge2.z);
      soa.nx_.push_back(normal.x);
      soa.ny_.push_back(normal.y);
      soa.nz_.push_back(normal.z);
      soa.min_x_.push_back(minimum.x);
      soa.min_y_.push_back(minimum.y);
      soa.min_z_.push_back(minimum.z);
      soa.max_x_.push_back(maximum.x);
      soa.max_y_.push_back(maximum.y);
      soa.max_z_.push_back(maximum.z);
      const auto vertices = canonical_vertices(triangle);
      Hash64 oriented_hash;
      oriented_hash.integer(triangle_key(
          vertices[0], subtract(vertices[1], vertices[0]),
          subtract(vertices[2], vertices[0])));
      oriented_hash.real(normal.x);
      oriented_hash.real(normal.y);
      oriented_hash.real(normal.z);
      canonical_triangle_hashes.push_back(oriented_hash.finish());
    }
    std::sort(canonical_triangle_hashes.begin(),
              canonical_triangle_hashes.end());
    Hash64 triangle_hash;
    triangle_hash.integer(static_cast<std::uint64_t>(triangle_count));
    for (const std::uint64_t value : canonical_triangle_hashes) {
      triangle_hash.integer(value);
    }
    soa.fingerprint_ = triangle_hash.finish();

    candidate.line_offsets_.assign(line_count + 1U, 0U);

    const Span<const double> first_global =
        geometry.axis(static_cast<CartesianAxis>(transverse[0])).centres();
    const Span<const double> second_global =
        geometry.axis(static_cast<CartesianAxis>(transverse[1])).centres();
    const std::size_t first_begin = static_cast<std::size_t>(
        component(patch.begin, transverse[0]));
    const std::size_t second_begin = static_cast<std::size_t>(
        component(patch.begin, transverse[1]));
    const double* first_centres = first_global.data + first_begin;
    const double* second_centres = second_global.data + second_begin;

    auto minimum_for = [&](int axis, std::size_t triangle) noexcept {
      return axis == 0 ? soa.min_x_[triangle]
                       : (axis == 1 ? soa.min_y_[triangle]
                                    : soa.min_z_[triangle]);
    };
    auto maximum_for = [&](int axis, std::size_t triangle) noexcept {
      return axis == 0 ? soa.max_x_[triangle]
                       : (axis == 1 ? soa.max_y_[triangle]
                                    : soa.max_z_[triangle]);
    };

    std::vector<std::size_t> bin_offsets(line_count + 1U, 0U);
    std::size_t reference_count = 0U;
    if (triangle_count != 0U) {
      for (std::size_t triangle = 0U; triangle < triangle_count; ++triangle) {
        const double first_min = minimum_for(transverse[0], triangle);
        const double first_max = maximum_for(transverse[0], triangle);
        const double second_min = minimum_for(transverse[1], triangle);
        const double second_max = maximum_for(transverse[1], triangle);
        const auto first_lo = std::lower_bound(first_centres,
                                               first_centres + first_count,
                                               first_min);
        const auto first_hi = std::upper_bound(first_centres,
                                               first_centres + first_count,
                                               first_max);
        const auto second_lo = std::lower_bound(second_centres,
                                                second_centres + second_count,
                                                second_min);
        const auto second_hi = std::upper_bound(second_centres,
                                                second_centres + second_count,
                                                second_max);
        const std::size_t lo0 = static_cast<std::size_t>(first_lo - first_centres);
        const std::size_t hi0 = static_cast<std::size_t>(first_hi - first_centres);
        const std::size_t lo1 =
            static_cast<std::size_t>(second_lo - second_centres);
        const std::size_t hi1 =
            static_cast<std::size_t>(second_hi - second_centres);
        std::size_t covered{};
        const std::size_t reference_limit = static_cast<std::size_t>(
            bin_reference_limit(budget));
        if (!checked_product(hi0 - lo0, hi1 - lo1, covered) ||
            reference_count > reference_limit ||
            covered > reference_limit - reference_count) {
          return invalid_plan(detail::stl_detail_bin_references);
        }
        reference_count += covered;
        for (std::size_t second = lo1; second < hi1; ++second) {
          for (std::size_t first = lo0; first < hi0; ++first) {
            ++bin_offsets[first + first_count * second + 1U];
          }
        }
      }
    }

    std::vector<std::size_t> bin_triangles;
    if (triangle_count != 0U) {
      for (std::size_t line = 0U; line < line_count; ++line) {
        bin_offsets[line + 1U] += bin_offsets[line];
      }
      bin_triangles.resize(reference_count);
      std::vector<std::size_t> cursor = bin_offsets;
      for (std::size_t triangle = 0U; triangle < triangle_count; ++triangle) {
        const auto first_lo = std::lower_bound(
            first_centres, first_centres + first_count,
            minimum_for(transverse[0], triangle));
        const auto first_hi = std::upper_bound(
            first_centres, first_centres + first_count,
            maximum_for(transverse[0], triangle));
        const auto second_lo = std::lower_bound(
            second_centres, second_centres + second_count,
            minimum_for(transverse[1], triangle));
        const auto second_hi = std::upper_bound(
            second_centres, second_centres + second_count,
            maximum_for(transverse[1], triangle));
        for (std::size_t second =
                 static_cast<std::size_t>(second_lo - second_centres);
             second < static_cast<std::size_t>(second_hi - second_centres);
             ++second) {
          for (std::size_t first =
                   static_cast<std::size_t>(first_lo - first_centres);
               first < static_cast<std::size_t>(first_hi - first_centres);
               ++first) {
            const std::size_t line = first + first_count * second;
            bin_triangles[cursor[line]++] = triangle;
          }
        }
      }
    }

    auto soa_real = [&](int vector_component, std::size_t triangle,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& z) noexcept {
      return vector_component == 0
                 ? x[triangle]
                 : (vector_component == 1 ? y[triangle] : z[triangle]);
    };
    struct WorkerOutput {
      std::vector<std::size_t> line_sizes;
      std::vector<double> span_bounds;
      Status status{};
    };
    const std::size_t worker_count = std::min<std::size_t>(
        line_count == 0U ? 1U : line_count,
        static_cast<std::size_t>(budget.worker_threads));
    std::vector<WorkerOutput> outputs(worker_count);
    const auto scan_chunk = [&](std::size_t worker) noexcept {
      WorkerOutput& output = outputs[worker];
      try {
        const detail::LineChunk chunk =
            detail::fixed_line_chunk(line_count, worker, worker_count);
        const std::size_t line_begin = chunk.begin;
        const std::size_t line_end = chunk.end;
        output.line_sizes.assign(line_end - line_begin, 0U);
        std::vector<Event> events;
        events.reserve(static_cast<std::size_t>(budget.max_events_per_line));
        const std::size_t chunk_candidate_references =
            bin_offsets[line_end] - bin_offsets[line_begin];
        output.span_bounds.reserve(std::min<std::size_t>(
            chunk_candidate_references,
            static_cast<std::size_t>(bin_reference_limit(budget))));
        for (std::size_t line = line_begin; line < line_end; ++line) {
        events.clear();
        const std::size_t first = line % first_count;
        const std::size_t second = line / first_count;
        Real3 origin{};
        if (transverse[0] == 0) {
          origin.x = first_centres[first];
        } else if (transverse[0] == 1) {
          origin.y = first_centres[first];
        } else {
          origin.z = first_centres[first];
        }
        if (transverse[1] == 0) {
          origin.x = second_centres[second];
        } else if (transverse[1] == 1) {
          origin.y = second_centres[second];
        } else {
          origin.z = second_centres[second];
        }

        const std::size_t candidate_begin = bin_offsets[line];
        const std::size_t candidate_end = bin_offsets[line + 1U];
        if (candidate_end - candidate_begin > budget.max_events_per_line) {
          output.status = invalid_plan(detail::stl_detail_budget);
          return;
        }
        for (std::size_t position = candidate_begin; position < candidate_end;
             ++position) {
          const std::size_t triangle = bin_triangles[position];
          const double coordinate0 = component(origin, transverse[0]);
          const double coordinate1 = component(origin, transverse[1]);
          if (coordinate0 < minimum_for(transverse[0], triangle) ||
              coordinate0 > maximum_for(transverse[0], triangle) ||
              coordinate1 < minimum_for(transverse[1], triangle) ||
              coordinate1 > maximum_for(transverse[1], triangle)) {
            continue;
          }
          const Real3 vertex{
              soa.ax_[triangle], soa.ay_[triangle], soa.az_[triangle]};
          const Real3 edge1{soa.e1x_[triangle], soa.e1y_[triangle],
                            soa.e1z_[triangle]};
          const Real3 edge2{soa.e2x_[triangle], soa.e2y_[triangle],
                            soa.e2z_[triangle]};
          Real3 direction{};
          if (scan == 0) {
            direction.x = 1.0;
          } else if (scan == 1) {
            direction.y = 1.0;
          } else {
            direction.z = 1.0;
          }
          const Real3 p = cross(direction, edge2);
          const double determinant = dot(edge1, p);
          const double determinant_tolerance =
              kRoundoffScale * std::numeric_limits<double>::epsilon() *
              std::sqrt(norm_squared(edge1) * norm_squared(edge2));
          if (std::abs(determinant) <= determinant_tolerance) {
            continue;
          }
          const Real3 displacement = subtract(origin, vertex);
          const double inverse_determinant = 1.0 / determinant;
          const double barycentric1 = dot(displacement, p) * inverse_determinant;
          const Real3 q = cross(displacement, edge1);
          const double barycentric2 = dot(direction, q) * inverse_determinant;
          constexpr double kBarycentricTolerance =
              kRoundoffScale * std::numeric_limits<double>::epsilon();
          if (barycentric1 < -kBarycentricTolerance ||
              barycentric2 < -kBarycentricTolerance ||
              barycentric1 + barycentric2 > 1.0 + kBarycentricTolerance) {
            continue;
          }
          const double intersection = dot(edge2, q) * inverse_determinant;
          if (!std::isfinite(intersection)) {
            output.status = invalid_plan(detail::stl_detail_triangle);
            return;
          }
          const double normal_component = soa_real(
              scan, triangle, soa.nx_, soa.ny_, soa.nz_);
          const TriangleInput original_triangle{
              vertex,
              {vertex.x + edge1.x, vertex.y + edge1.y, vertex.z + edge1.z},
              {vertex.x + edge2.x, vertex.y + edge2.y, vertex.z + edge2.z}};
          const auto canonical = canonical_vertices(original_triangle);
          if (events.size() >= budget.max_events_per_line) {
            output.status = invalid_plan(detail::stl_detail_budget);
            return;
          }
          events.push_back({intersection,
                            static_cast<std::int8_t>(normal_component > 0.0
                                                         ? 1
                                                         : -1),
                            triangle_key(
                                canonical[0],
                                subtract(canonical[1], canonical[0]),
                                subtract(canonical[2], canonical[0]))});
        }

        std::sort(events.begin(), events.end(), event_less);
        std::size_t write = 0U;
        for (std::size_t read = 0U; read < events.size();) {
          const double anchor = events[read].coordinate;
          bool kept_negative = false;
          bool kept_positive = false;
          std::size_t group_end = read;
          while (group_end < events.size() &&
                 same_coordinate(anchor, events[group_end].coordinate)) {
            const bool positive = events[group_end].normal_sign > 0;
            bool& kept = positive ? kept_positive : kept_negative;
            if (!kept) {
              events[write++] = events[group_end];
              kept = true;
            }
            ++group_end;
          }
          read = group_end;
        }
        events.resize(write);
        if ((events.size() & 1U) != 0U) {
          output.status = invalid_plan(detail::stl_detail_open_surface);
          return;
        }
        for (std::size_t index = 0U; index < events.size(); index += 2U) {
          const std::uint64_t span_limit = bin_reference_limit(budget);
          if (span_limit < 2U || output.span_bounds.size() > span_limit - 2U) {
            output.status = invalid_plan(detail::stl_detail_budget);
            return;
          }
          output.span_bounds.push_back(events[index].coordinate);
          output.span_bounds.push_back(events[index + 1U].coordinate);
        }
        output.line_sizes[line - line_begin] = events.size();
        }
      } catch (const std::bad_alloc&) {
        output.status = allocation_failure();
      } catch (...) {
        output.status = invalid_plan(detail::stl_detail_triangle);
      }
    };
    std::vector<std::thread> workers;
    workers.reserve(worker_count > 0U ? worker_count - 1U : 0U);
    Status launch_status{};
    try {
      for (std::size_t worker = 1U; worker < worker_count; ++worker) {
        if (inject_thread_launch_failure()) {
          throw std::bad_alloc{};
        }
        workers.emplace_back(scan_chunk, worker);
      }
    } catch (...) {
      launch_status = allocation_failure();
    }
    if (launch_status) {
      scan_chunk(0U);
    }
    for (std::thread& worker : workers) {
      worker.join();
    }
    if (!launch_status) {
      return launch_status;
    }
    for (const WorkerOutput& output : outputs) {
      if (!output.status) {
        return output.status;
      }
    }
    std::size_t total_span_bounds = 0U;
    const std::size_t span_limit =
        static_cast<std::size_t>(bin_reference_limit(budget));
    for (const WorkerOutput& output : outputs) {
      if (total_span_bounds > span_limit ||
          output.span_bounds.size() > span_limit - total_span_bounds) {
        return invalid_plan(detail::stl_detail_budget);
      }
      total_span_bounds += output.span_bounds.size();
    }
    candidate.span_bounds_.reserve(total_span_bounds);
    std::size_t line = 0U;
    for (WorkerOutput& output : outputs) {
      std::size_t cumulative = candidate.span_bounds_.size();
      for (const std::size_t line_size : output.line_sizes) {
        cumulative += line_size;
        candidate.line_offsets_[++line] = cumulative;
      }
      candidate.span_bounds_.insert(candidate.span_bounds_.end(),
                                    output.span_bounds.begin(),
                                    output.span_bounds.end());
    }

    Hash64 plan_hash;
    plan_hash.integer(geometry.fingerprint());
    plan_hash.integer(static_cast<std::uint8_t>(scan));
    plan_hash.integer(patch.begin.x);
    plan_hash.integer(patch.begin.y);
    plan_hash.integer(patch.begin.z);
    plan_hash.integer(patch.cells.x);
    plan_hash.integer(patch.cells.y);
    plan_hash.integer(patch.cells.z);
    plan_hash.integer(soa.fingerprint());
    plan_hash.integer(static_cast<std::uint64_t>(candidate.line_offsets_.size()));
    for (const std::size_t offset : candidate.line_offsets_) {
      plan_hash.integer(static_cast<std::uint64_t>(offset));
    }
    plan_hash.integer(static_cast<std::uint64_t>(candidate.span_bounds_.size()));
    for (const double bound : candidate.span_bounds_) {
      plan_hash.real(bound);
    }
    candidate.fingerprint_ = plan_hash.finish();
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return allocation_failure();
  } catch (...) {
    return invalid_plan(detail::stl_detail_triangle);
  }
}

Status StlScanPlan::classify(const CartesianGeometryPlan& geometry,
                             const MeshPatch& patch,
                             Span<std::uint8_t> region) const noexcept {
  const int scan = axis_index(scan_axis_);
  if (scan < 0 || fingerprint_ == 0U || geometry.fingerprint() == 0U ||
      geometry.fingerprint() != geometry_fingerprint_ ||
      patch.begin.x != patch_begin_.x || patch.begin.y != patch_begin_.y ||
      patch.begin.z != patch_begin_.z || patch.cells.x != patch_cells_.x ||
      patch.cells.y != patch_cells_.y || patch.cells.z != patch_cells_.z) {
    return invalid_plan(detail::stl_detail_geometry);
  }
  if (patch.cells.x <= 0 || patch.cells.y <= 0 || patch.cells.z <= 0) {
    return invalid_plan(detail::stl_detail_patch);
  }
  std::size_t xy{};
  std::size_t cell_count{};
  if (!checked_product(static_cast<std::size_t>(patch.cells.x),
                       static_cast<std::size_t>(patch.cells.y), xy) ||
      !checked_product(xy, static_cast<std::size_t>(patch.cells.z),
                       cell_count) ||
      region.size != cell_count || (cell_count != 0U && region.data == nullptr)) {
    return invalid_plan(detail::stl_detail_region);
  }
  const auto transverse = transverse_axes(scan);
  std::size_t line_count{};
  if (!checked_product(
          static_cast<std::size_t>(component(patch.cells, transverse[0])),
          static_cast<std::size_t>(component(patch.cells, transverse[1])),
          line_count) ||
      line_offsets_.size() != line_count + 1U ||
      line_offsets_.back() != span_bounds_.size()) {
    return invalid_plan(detail::stl_detail_region);
  }

  std::fill(region.data, region.data + region.size,
            static_cast<std::uint8_t>(RegionFlag::fluid));
  const Span<const double> scan_centres_global =
      geometry.axis(scan_axis_).centres();
  const std::size_t scan_begin =
      static_cast<std::size_t>(component(patch.begin, scan));
  const std::size_t scan_count =
      static_cast<std::size_t>(component(patch.cells, scan));
  if (scan_begin > scan_centres_global.size ||
      scan_count > scan_centres_global.size - scan_begin) {
    return invalid_plan(detail::stl_detail_geometry);
  }
  const double* scan_centres = scan_centres_global.data + scan_begin;
  const std::size_t first_count = static_cast<std::size_t>(
      component(patch.cells, transverse[0]));
  for (std::size_t line = 0U; line < line_count; ++line) {
    const std::size_t first = line % first_count;
    const std::size_t second = line / first_count;
    std::size_t span = line_offsets_[line];
    const std::size_t span_end = line_offsets_[line + 1U];
    if (((span_end - span) & 1U) != 0U || span_end > span_bounds_.size()) {
      return invalid_plan(detail::stl_detail_region);
    }
    for (std::size_t along = 0U; along < scan_count; ++along) {
      const double coordinate_value = scan_centres[along];
      while (span < span_end && coordinate_value > span_bounds_[span + 1U]) {
        span += 2U;
      }
      if (span >= span_end || coordinate_value < span_bounds_[span]) {
        continue;
      }
      std::size_t x{}, y{}, z{};
      if (scan == 0) {
        x = along;
        y = first;
        z = second;
      } else if (scan == 1) {
        x = first;
        y = along;
        z = second;
      } else {
        x = first;
        y = second;
        z = along;
      }
      const std::size_t flat =
          x + static_cast<std::size_t>(patch.cells.x) *
                  (y + static_cast<std::size_t>(patch.cells.y) * z);
      region.data[flat] = static_cast<std::uint8_t>(RegionFlag::solid);
    }
  }
  return {};
}

Status StlScanCompiler::compile(
    MPI_Comm communicator, const std::filesystem::path& case_root,
    const std::optional<std::filesystem::path>& stl_file,
    const CartesianGeometryPlan& geometry, const MeshPatch& patch,
    CartesianAxis scan_axis, StlScanBudget budget,
    StlScanPlan& out) noexcept {
  if (communicator == MPI_COMM_NULL) {
    return mpi_failure();
  }
  int rank = -1;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS) {
    return mpi_failure();
  }

  const int scan = axis_index(scan_axis);
  Status local_input{};
  std::size_t local_line_count = 0U;
  if (scan < 0 || !valid_budget(budget)) {
    local_input = invalid_plan(detail::stl_detail_budget);
  } else {
    const auto transverse = transverse_axes(scan);
    const std::int32_t first = component(patch.cells, transverse[0]);
    const std::int32_t second = component(patch.cells, transverse[1]);
    if (first <= 0 || second <= 0 ||
        !checked_product(static_cast<std::size_t>(first),
                         static_cast<std::size_t>(second), local_line_count)) {
      local_input = invalid_plan(detail::stl_detail_patch);
    } else {
      local_input = initial_memory_gate(0U, local_line_count, 0U, budget);
    }
  }
  const Status input_consensus = collective_status(communicator, rank,
                                                   local_input);
  if (!input_consensus) {
    return input_consensus;
  }

  std::vector<TriangleInput> triangles;
  std::uint64_t source_bytes = 0U;
  Status root_status{};
  try {
    if (rank == 0 && stl_file.has_value()) {
      std::vector<std::uint8_t> bytes;
      std::size_t triangle_upper_bound = 0U;
      root_status = read_stl_file(case_root, *stl_file, local_line_count,
                                  budget, bytes, triangle_upper_bound);
      if (root_status) {
        source_bytes = static_cast<std::uint64_t>(bytes.size());
        triangles.reserve(triangle_upper_bound);
        root_status = detail::parse_stl_bytes(
            {bytes.data(), bytes.size()}, triangles);
      }
    }
  } catch (const std::bad_alloc&) {
    root_status = allocation_failure();
  } catch (...) {
    root_status = io_failure(detail::stl_detail_read);
  }

  std::array<std::uint32_t, 2U> root_wire{};
  if (rank == 0) {
    root_wire = {static_cast<std::uint32_t>(root_status.code),
                 root_status.detail};
  }
  if (MPI_Bcast(root_wire.data(), static_cast<int>(root_wire.size()),
                MPI_UINT32_T, 0, communicator) != MPI_SUCCESS) {
    return mpi_failure();
  }
  root_status = {static_cast<StatusCode>(root_wire[0]), root_wire[1]};
  if (!root_status) {
    return root_status;
  }

  std::array<std::uint64_t, 2U> input_header{
      rank == 0 ? static_cast<std::uint64_t>(triangles.size()) : 0U,
      rank == 0 ? source_bytes : 0U};
  if (MPI_Bcast(input_header.data(), static_cast<int>(input_header.size()),
                MPI_UINT64_T, 0, communicator) !=
      MPI_SUCCESS) {
    return mpi_failure();
  }
  const std::uint64_t triangle_count = input_header[0];
  source_bytes = input_header[1];
  Status resize_status{};
  if (triangle_count > detail::kMaxStlTriangles ||
      triangle_count > std::numeric_limits<std::size_t>::max()) {
    resize_status = invalid_plan(detail::stl_detail_triangle);
  } else {
    if (scan < 0) {
      resize_status = invalid_plan(detail::stl_detail_geometry);
    } else {
      const auto transverse = transverse_axes(scan);
      std::size_t line_count{};
      if (!checked_product(
              static_cast<std::size_t>(component(patch.cells, transverse[0])),
              static_cast<std::size_t>(component(patch.cells, transverse[1])),
              line_count)) {
        resize_status = invalid_plan(detail::stl_detail_patch);
      } else {
        resize_status = initial_memory_gate(
            static_cast<std::size_t>(triangle_count), line_count,
            rank == 0 ? source_bytes : 0U, budget);
      }
    }
  }
  if (resize_status && rank != 0) {
    try {
      triangles.resize(static_cast<std::size_t>(triangle_count));
    } catch (const std::bad_alloc&) {
      resize_status = allocation_failure();
    } catch (...) {
      resize_status = allocation_failure();
    }
  }
  const Status resize_collective =
      collective_status(communicator, rank, resize_status);
  if (!resize_collective) {
    return resize_collective;
  }
  std::vector<double> coordinate_wire;
  Status wire_status{};
  try {
    std::size_t coordinate_count{};
    if (!checked_product(triangles.size(), 9U, coordinate_count)) {
      wire_status = invalid_plan(detail::stl_detail_triangle);
    } else {
      coordinate_wire.resize(coordinate_count);
      if (rank == 0) {
        for (std::size_t index = 0U; index < triangles.size(); ++index) {
          const TriangleInput& triangle = triangles[index];
          const std::array<double, 9U> coordinates{
              triangle.a.x, triangle.a.y, triangle.a.z,
              triangle.b.x, triangle.b.y, triangle.b.z,
              triangle.c.x, triangle.c.y, triangle.c.z};
          std::copy(coordinates.begin(), coordinates.end(),
                    coordinate_wire.data() + index * coordinates.size());
        }
      }
    }
  } catch (...) {
    wire_status = allocation_failure();
  }
  const Status wire_collective =
      collective_status(communicator, rank, wire_status);
  if (!wire_collective) {
    return wire_collective;
  }
  const Status broadcast_status = broadcast_doubles(
      communicator, 0, coordinate_wire.data(), coordinate_wire.size());
  if (!broadcast_status) {
    return broadcast_status;
  }
  if (rank != 0) {
    for (std::size_t index = 0U; index < triangles.size(); ++index) {
      const double* coordinates = coordinate_wire.data() + index * 9U;
      triangles[index] = {{coordinates[0], coordinates[1], coordinates[2]},
                          {coordinates[3], coordinates[4], coordinates[5]},
                          {coordinates[6], coordinates[7], coordinates[8]}};
    }
  }

  StlScanPlan candidate;
  const Status local = compile_triangles(geometry, patch,
                                         {triangles.data(), triangles.size()},
                                         scan_axis, budget, candidate);
  const Status compiled = collective_status(communicator, rank, local);
  if (!compiled) {
    return compiled;
  }
  out = std::move(candidate);
  return {};
}

}  // namespace hundun::v04
