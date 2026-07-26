// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/performance_correctness.hpp"
#include "yyjson.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

struct Sha256 final {
  std::array<std::uint32_t, 8> state{
      UINT32_C(0x6a09e667), UINT32_C(0xbb67ae85),
      UINT32_C(0x3c6ef372), UINT32_C(0xa54ff53a),
      UINT32_C(0x510e527f), UINT32_C(0x9b05688c),
      UINT32_C(0x1f83d9ab), UINT32_C(0x5be0cd19)};
  std::array<std::uint8_t, 64> block{};
  std::uint64_t bytes{};
  std::size_t used{};
};

constexpr std::array<std::uint32_t, 64> kSha256{
    UINT32_C(0x428a2f98), UINT32_C(0x71374491), UINT32_C(0xb5c0fbcf),
    UINT32_C(0xe9b5dba5), UINT32_C(0x3956c25b), UINT32_C(0x59f111f1),
    UINT32_C(0x923f82a4), UINT32_C(0xab1c5ed5), UINT32_C(0xd807aa98),
    UINT32_C(0x12835b01), UINT32_C(0x243185be), UINT32_C(0x550c7dc3),
    UINT32_C(0x72be5d74), UINT32_C(0x80deb1fe), UINT32_C(0x9bdc06a7),
    UINT32_C(0xc19bf174), UINT32_C(0xe49b69c1), UINT32_C(0xefbe4786),
    UINT32_C(0x0fc19dc6), UINT32_C(0x240ca1cc), UINT32_C(0x2de92c6f),
    UINT32_C(0x4a7484aa), UINT32_C(0x5cb0a9dc), UINT32_C(0x76f988da),
    UINT32_C(0x983e5152), UINT32_C(0xa831c66d), UINT32_C(0xb00327c8),
    UINT32_C(0xbf597fc7), UINT32_C(0xc6e00bf3), UINT32_C(0xd5a79147),
    UINT32_C(0x06ca6351), UINT32_C(0x14292967), UINT32_C(0x27b70a85),
    UINT32_C(0x2e1b2138), UINT32_C(0x4d2c6dfc), UINT32_C(0x53380d13),
    UINT32_C(0x650a7354), UINT32_C(0x766a0abb), UINT32_C(0x81c2c92e),
    UINT32_C(0x92722c85), UINT32_C(0xa2bfe8a1), UINT32_C(0xa81a664b),
    UINT32_C(0xc24b8b70), UINT32_C(0xc76c51a3), UINT32_C(0xd192e819),
    UINT32_C(0xd6990624), UINT32_C(0xf40e3585), UINT32_C(0x106aa070),
    UINT32_C(0x19a4c116), UINT32_C(0x1e376c08), UINT32_C(0x2748774c),
    UINT32_C(0x34b0bcb5), UINT32_C(0x391c0cb3), UINT32_C(0x4ed8aa4a),
    UINT32_C(0x5b9cca4f), UINT32_C(0x682e6ff3), UINT32_C(0x748f82ee),
    UINT32_C(0x78a5636f), UINT32_C(0x84c87814), UINT32_C(0x8cc70208),
    UINT32_C(0x90befffa), UINT32_C(0xa4506ceb), UINT32_C(0xbef9a3f7),
    UINT32_C(0xc67178f2)};

std::uint32_t rotate_right(std::uint32_t value, unsigned amount) noexcept {
  return (value >> amount) | (value << (32U - amount));
}

void sha256_transform(Sha256& hash) noexcept {
  std::array<std::uint32_t, 64> words{};
  for (std::size_t index = 0; index < 16U; ++index) {
    const auto offset = index * 4U;
    words[index] =
        (static_cast<std::uint32_t>(hash.block[offset]) << 24U) |
        (static_cast<std::uint32_t>(hash.block[offset + 1U]) << 16U) |
        (static_cast<std::uint32_t>(hash.block[offset + 2U]) << 8U) |
        static_cast<std::uint32_t>(hash.block[offset + 3U]);
  }
  for (std::size_t index = 16U; index < words.size(); ++index) {
    const auto x = words[index - 15U];
    const auto y = words[index - 2U];
    const auto s0 =
        rotate_right(x, 7U) ^ rotate_right(x, 18U) ^ (x >> 3U);
    const auto s1 =
        rotate_right(y, 17U) ^ rotate_right(y, 19U) ^ (y >> 10U);
    words[index] =
        words[index - 16U] + s0 + words[index - 7U] + s1;
  }
  auto a = hash.state[0U];
  auto b = hash.state[1U];
  auto c = hash.state[2U];
  auto d = hash.state[3U];
  auto e = hash.state[4U];
  auto f = hash.state[5U];
  auto g = hash.state[6U];
  auto h = hash.state[7U];
  for (std::size_t index = 0; index < words.size(); ++index) {
    const auto sum1 = rotate_right(e, 6U) ^ rotate_right(e, 11U) ^
                      rotate_right(e, 25U);
    const auto choice = (e & f) ^ (~e & g);
    const auto first = h + sum1 + choice + kSha256[index] + words[index];
    const auto sum0 = rotate_right(a, 2U) ^ rotate_right(a, 13U) ^
                      rotate_right(a, 22U);
    const auto majority = (a & b) ^ (a & c) ^ (b & c);
    const auto second = sum0 + majority;
    h = g;
    g = f;
    f = e;
    e = d + first;
    d = c;
    c = b;
    b = a;
    a = first + second;
  }
  hash.state[0U] += a;
  hash.state[1U] += b;
  hash.state[2U] += c;
  hash.state[3U] += d;
  hash.state[4U] += e;
  hash.state[5U] += f;
  hash.state[6U] += g;
  hash.state[7U] += h;
}

void sha256_update(Sha256& hash, std::string_view bytes) noexcept {
  for (const char character : bytes) {
    const auto byte = static_cast<unsigned char>(character);
    hash.block[hash.used++] = byte;
    ++hash.bytes;
    if (hash.used == hash.block.size()) {
      sha256_transform(hash);
      hash.used = 0U;
    }
  }
}

std::string sha256(std::string_view bytes) {
  Sha256 hash;
  sha256_update(hash, bytes);
  const std::uint64_t bit_count = hash.bytes * 8U;
  hash.block[hash.used++] = UINT8_C(0x80);
  if (hash.used > 56U) {
    while (hash.used < hash.block.size())
      hash.block[hash.used++] = 0U;
    sha256_transform(hash);
    hash.used = 0U;
  }
  while (hash.used < 56U)
    hash.block[hash.used++] = 0U;
  for (int shift = 56; shift >= 0; shift -= 8)
    hash.block[hash.used++] =
        static_cast<std::uint8_t>(bit_count >> shift);
  sha256_transform(hash);
  std::ostringstream output;
  output << std::hex << std::setfill('0');
  for (const auto value : hash.state)
    output << std::setw(8) << value;
  return output.str();
}

std::string read_file(const std::filesystem::path& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream)
    throw std::runtime_error("unable to read evidence input file");
  return {std::istreambuf_iterator<char>(stream),
          std::istreambuf_iterator<char>()};
}

yyjson_val* required(yyjson_val* object, const char* key) {
  auto* value = yyjson_obj_get(object, key);
  if (value == nullptr)
    throw std::runtime_error(std::string("missing JSON member: ") + key);
  return value;
}

std::uint64_t u64(yyjson_val* object, const char* key) {
  auto* value = required(object, key);
  if (!yyjson_is_uint(value))
    throw std::runtime_error(std::string("JSON member is not u64: ") + key);
  return yyjson_get_uint(value);
}

double real(yyjson_val* object, const char* key) {
  auto* value = required(object, key);
  if (!yyjson_is_num(value))
    throw std::runtime_error(std::string("JSON member is not real: ") + key);
  const double result = yyjson_get_real(value);
  if (!std::isfinite(result) || result < 0.0)
    throw std::runtime_error(std::string("invalid real JSON member: ") + key);
  return result;
}

std::string string_value(yyjson_val* object, const char* key) {
  auto* value = required(object, key);
  if (!yyjson_is_str(value))
    throw std::runtime_error(std::string("JSON member is not string: ") + key);
  return yyjson_get_str(value);
}

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

struct Halo final {
  std::array<std::uint64_t, 9> integer{};
  double wait{};
};

Halo halo(yyjson_val* value) {
  Halo result;
  constexpr std::array<const char*, 9> names{
      "completed_exchanges", "begin_calls", "wait_calls",
      "send_payload_bytes", "receive_payload_bytes", "pack_bytes",
      "unpack_bytes", "send_messages", "receive_messages"};
  for (std::size_t index = 0; index < names.size(); ++index)
    result.integer[index] = u64(value, names[index]);
  result.wait = real(value, "successful_wait_seconds");
  return result;
}

Halo add(Halo left, const Halo& right) {
  for (std::size_t index = 0; index < left.integer.size(); ++index) {
    if (right.integer[index] >
        std::numeric_limits<std::uint64_t>::max() - left.integer[index])
      throw std::runtime_error("Halo evidence integer overflow");
    left.integer[index] += right.integer[index];
  }
  left.wait += right.wait;
  if (!std::isfinite(left.wait))
    throw std::runtime_error("Halo evidence wait overflow");
  return left;
}

std::uint64_t checked_add(std::uint64_t left, std::uint64_t right) {
  if (right > std::numeric_limits<std::uint64_t>::max() - left)
    throw std::runtime_error("performance evidence integer overflow");
  return left + right;
}

std::uint64_t checked_multiply(std::uint64_t left, std::uint64_t right) {
  if (left != 0U &&
      right > std::numeric_limits<std::uint64_t>::max() / left)
    throw std::runtime_error("performance evidence coverage overflow");
  return left * right;
}

bool equal(const Halo& left, const Halo& right) noexcept {
  return left.integer == right.integer && bits(left.wait) == bits(right.wait);
}

struct Entry final {
  std::uint64_t repetition{};
  std::uint64_t rank{};
  Halo runtime;
  Halo pressure;
};

std::array<std::uint64_t, 3> u64_triplet(yyjson_val* object,
                                         const char* key) {
  auto* value = required(object, key);
  if (!yyjson_is_arr(value) || yyjson_arr_size(value) != 3U)
    throw std::runtime_error(std::string("JSON triplet differs: ") + key);
  std::array<std::uint64_t, 3> result{};
  for (std::size_t index = 0; index < result.size(); ++index) {
    auto* item = yyjson_arr_get(value, index);
    if (!yyjson_is_uint(item) || yyjson_get_uint(item) == 0U)
      throw std::runtime_error(std::string("JSON triplet is invalid: ") + key);
    result[index] = yyjson_get_uint(item);
  }
  return result;
}

void validate_fixture_halo_oracle(
    const std::vector<Entry>& entries, yyjson_val* artifact_root) {
  constexpr std::uint64_t runtime_calls_per_step = 32U;
  constexpr std::uint64_t runtime_components_per_step = 98U;
  constexpr std::uint64_t pressure_applies_per_step = 6U;
  constexpr std::uint64_t halo_depth = 2U;
  constexpr std::uint64_t runtime_region_count = 26U;
  const auto steps =
      u64(required(artifact_root, "measurement"), "measured_steps");
  const auto local =
      u64_triplet(artifact_root, "per_rank_owned_cell_extents");
  const auto grid = u64_triplet(artifact_root, "process_grid");
  const auto owned_plane =
      checked_multiply(local[0U], local[1U]);
  const auto owned = checked_multiply(owned_plane, local[2U]);
  std::uint64_t extended = 1U;
  for (const auto extent : local)
    extended = checked_multiply(
        extended,
        checked_add(extent, checked_multiply(halo_depth, 2U)));
  if (extended < owned)
    throw std::runtime_error("runtime fixture region underflow");
  const auto runtime_region_values = extended - owned;
  const auto runtime_calls =
      checked_multiply(steps, runtime_calls_per_step);
  const auto runtime_components =
      checked_multiply(steps, runtime_components_per_step);
  const auto runtime_bytes = checked_multiply(
      checked_multiply(runtime_region_values, runtime_components),
      sizeof(double));
  const auto runtime_messages =
      checked_multiply(runtime_region_count, runtime_calls);

  std::uint64_t pressure_values = 0U;
  std::uint64_t pressure_neighbors = 0U;
  for (std::size_t axis = 0U; axis < grid.size(); ++axis) {
    if (grid[axis] == 1U)
      continue;
    std::uint64_t face_values = 1U;
    for (std::size_t other = 0U; other < local.size(); ++other)
      if (other != axis)
        face_values = checked_multiply(face_values, local[other]);
    pressure_values =
        checked_add(pressure_values,
                    checked_multiply(face_values, 2U));
    pressure_neighbors =
        checked_add(pressure_neighbors, grid[axis] == 2U ? 1U : 2U);
  }
  const auto pressure_calls =
      checked_multiply(steps, pressure_applies_per_step);
  const auto pressure_bytes = checked_multiply(
      checked_multiply(pressure_values, pressure_calls), sizeof(double));
  const auto pressure_messages =
      checked_multiply(pressure_neighbors, pressure_calls);

  const auto exact = [](const Halo& value, std::uint64_t calls,
                        std::uint64_t bytes,
                        std::uint64_t messages) noexcept {
    return value.integer ==
           std::array<std::uint64_t, 9>{
               calls, calls, calls, bytes, bytes, bytes, bytes,
               messages, messages};
  };
  for (const auto& entry : entries) {
    if (!exact(entry.runtime, runtime_calls, runtime_bytes,
               runtime_messages) ||
        !exact(entry.pressure, pressure_calls, pressure_bytes,
               pressure_messages))
      throw std::runtime_error(
          "Halo evidence differs from the independent fixture oracle");
  }
}

void validate_metric(yyjson_val* metric, bool applicable, double expected,
                     const char* unit) {
  if (string_value(metric, "unit") != unit)
    throw std::runtime_error("performance evidence metric unit differs");
  const std::string expected_status =
      applicable ? "available" : "not_applicable";
  if (string_value(metric, "status") != expected_status)
    throw std::runtime_error("performance evidence metric status differs");
  auto* value = required(metric, "value");
  if (!applicable) {
    if (!yyjson_is_null(value))
      throw std::runtime_error("not-applicable metric has a value");
  } else if (!yyjson_is_num(value) ||
             bits(yyjson_get_real(value)) != bits(expected)) {
    throw std::runtime_error("performance evidence metric value differs");
  }
}

void validate_input(std::string_view input, std::string_view artifact) {
  yyjson_doc* input_doc = yyjson_read(input.data(), input.size(), 0);
  yyjson_doc* artifact_doc =
      yyjson_read(artifact.data(), artifact.size(), 0);
  if (input_doc == nullptr || artifact_doc == nullptr) {
    if (input_doc != nullptr)
      yyjson_doc_free(input_doc);
    if (artifact_doc != nullptr)
      yyjson_doc_free(artifact_doc);
    throw std::runtime_error("invalid performance evidence JSON");
  }
  auto free_documents = [&] {
    yyjson_doc_free(input_doc);
    yyjson_doc_free(artifact_doc);
  };
  try {
    auto* root = yyjson_doc_get_root(input_doc);
    auto* artifact_root = yyjson_doc_get_root(artifact_doc);
    if (u64(root, "schema_version") != 1U)
      throw std::runtime_error("unsupported evidence input schema");
    if (u64(artifact_root, "schema_version") != 1U)
      throw std::runtime_error("unsupported performance artifact schema");
    auto* correctness = required(artifact_root, "correctness");
    auto* passed = required(correctness, "passed");
    if (!yyjson_is_bool(passed) || !yyjson_get_bool(passed))
      throw std::runtime_error("performance artifact correctness failed");
    const auto correctness_record =
        hundun::diagnostics::parse_performance_correctness(
            string_value(correctness, "summary"));
    if (!correctness_record.passed)
      throw std::runtime_error(
          "performance correctness summary did not pass");
    auto* measurement = required(artifact_root, "measurement");
    hundun::diagnostics::validate_performance_correctness_coverage(
        correctness_record, u64(measurement, "warmup_steps"),
        u64(measurement, "measured_steps"),
        u64(measurement, "repetitions"));
    auto* entries_value = required(root, "entries");
    if (!yyjson_is_arr(entries_value) || yyjson_arr_size(entries_value) == 0U)
      throw std::runtime_error("performance evidence entries are empty");
    std::vector<Entry> entries;
    const std::size_t count = yyjson_arr_size(entries_value);
    entries.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
      auto* value = yyjson_arr_get(entries_value, index);
      Entry entry;
      entry.repetition = u64(value, "repetition");
      entry.rank = u64(value, "relative_rank");
      entry.runtime = halo(required(value, "runtime"));
      entry.pressure = halo(required(value, "pressure"));
      const auto combined = halo(required(value, "combined"));
      if (!equal(combined, add(entry.runtime, entry.pressure)))
        throw std::runtime_error("combined Halo evidence differs");
      entries.push_back(entry);
    }
    std::uint64_t inferred_ranks = 0U;
    std::uint64_t repetitions = 0U;
    for (const auto& entry : entries) {
      inferred_ranks =
          std::max(inferred_ranks, checked_add(entry.rank, 1U));
      repetitions =
          std::max(repetitions, checked_add(entry.repetition, 1U));
    }
    if (inferred_ranks == 0U || repetitions == 0U)
      throw std::runtime_error("performance evidence dimensions are zero");
    if (entries.size() !=
        static_cast<std::size_t>(
            checked_multiply(inferred_ranks, repetitions)))
      throw std::runtime_error("performance evidence coverage differs");
    if (u64(artifact_root, "ranks") != inferred_ranks ||
        u64(required(artifact_root, "measurement"), "repetitions") !=
            repetitions)
      throw std::runtime_error(
          "artifact rank/repetition binding differs from evidence");
    for (std::size_t index = 0; index < entries.size(); ++index) {
      if (entries[index].repetition != index / inferred_ranks ||
          entries[index].rank != index % inferred_ranks)
        throw std::runtime_error("performance evidence order differs");
    }
    validate_fixture_halo_oracle(entries, artifact_root);

    Halo total;
    for (const auto& entry : entries)
      total = add(total, add(entry.runtime, entry.pressure));
    auto* counters = required(artifact_root, "exact_counters");
    auto* payload = required(counters, "halo_payload_bytes");
    auto* messages = required(counters, "halo_messages");
    if (u64(payload, "send") != total.integer[3U] ||
        u64(payload, "receive") != total.integer[4U] ||
        u64(payload, "pack") != total.integer[5U] ||
        u64(payload, "unpack") != total.integer[6U] ||
        u64(messages, "send") != total.integer[7U] ||
        u64(messages, "receive") != total.integer[8U])
      throw std::runtime_error(
          "artifact Halo maps differ from independent evidence sum");

    auto* summaries = required(root, "halo_summaries");
    if (!yyjson_is_arr(summaries) ||
        yyjson_arr_size(summaries) !=
            checked_multiply(repetitions, 3U))
      throw std::runtime_error("Halo summary coverage differs");
    constexpr std::array<const char*, 3> names{
        "runtime", "pressure", "combined"};
    for (std::uint64_t repetition = 0U; repetition < repetitions;
         ++repetition) {
      for (std::size_t family = 0U; family < names.size(); ++family) {
        std::uint64_t receive = 0U;
        double critical = 0.0;
        for (const auto& entry : entries) {
          if (entry.repetition != repetition)
            continue;
          const auto value =
              family == 0U
                  ? entry.runtime
                  : (family == 1U ? entry.pressure
                                  : add(entry.runtime, entry.pressure));
          receive = checked_add(receive, value.integer[4U]);
          critical = std::max(critical, value.wait);
        }
        const bool applicable = receive != 0U;
        if (applicable && critical <= 0.0)
          throw std::runtime_error(
              "nonlocal Halo evidence has no critical wait");
        auto* summary = yyjson_arr_get(
            summaries,
            static_cast<std::size_t>(repetition) * 3U + family);
        if (u64(summary, "repetition") != repetition ||
            string_value(summary, "family") != names[family])
          throw std::runtime_error("Halo summary key differs");
        validate_metric(required(summary, "critical_wait_seconds"),
                        applicable, critical, "s");
        validate_metric(
            required(summary,
                     "effective_bandwidth_bytes_per_second"),
            applicable,
            applicable ? static_cast<double>(receive) / critical : 0.0,
            "byte/s");
      }
    }

    auto* io = required(root, "io");
    auto* logical_io = required(counters, "logical_io_bytes");
    for (const char* name : {"checkpoint", "diagnostics"}) {
      auto* item = required(io, name);
      const auto logical = u64(item, "logical_bytes");
      const auto actual = u64(item, "actual_bytes");
      if (logical != u64(logical_io, name))
        throw std::runtime_error("artifact logical I/O differs");
      if ((logical == 0U) != (actual == 0U))
        throw std::runtime_error(
            "logical and actual I/O applicability differs");
      const bool applicable = logical != 0U;
      auto* duration_metric = required(item, "duration_seconds");
      const double duration =
          applicable ? yyjson_get_real(required(duration_metric, "value"))
                     : 0.0;
      if (applicable && (!std::isfinite(duration) || duration <= 0.0))
        throw std::runtime_error("I/O duration is invalid");
      validate_metric(duration_metric, applicable, duration, "s");
      validate_metric(
          required(item, "throughput_bytes_per_second"), applicable,
          applicable ? static_cast<double>(actual) / duration : 0.0,
          "byte/s");
    }
  } catch (...) {
    free_documents();
    throw;
  }
  free_documents();
}

bool hex64(std::string_view value) noexcept {
  if (value.size() != 64U)
    return false;
  for (const char character : value)
    if (!((character >= '0' && character <= '9') ||
          (character >= 'a' && character <= 'f')))
      return false;
  return true;
}

std::string escape_json(std::string_view value) {
  std::string result;
  for (const char raw_character : value) {
    const auto character =
        static_cast<unsigned char>(raw_character);
    switch (character) {
      case '"':
        result += "\\\"";
        break;
      case '\\':
        result += "\\\\";
        break;
      case '\n':
        result += "\\n";
        break;
      case '\r':
        result += "\\r";
        break;
      case '\t':
        result += "\\t";
        break;
      default:
        if (character < 0x20U)
          throw std::runtime_error("unsupported control character");
        result.push_back(static_cast<char>(character));
    }
  }
  return result;
}

std::vector<std::string> split_tabs(std::string_view line) {
  std::vector<std::string> result;
  std::size_t begin = 0U;
  for (;;) {
    const auto end = line.find('\t', begin);
    result.emplace_back(line.substr(
        begin, end == std::string_view::npos ? line.size() - begin
                                             : end - begin));
    if (end == std::string_view::npos)
      return result;
    begin = end + 1U;
  }
}

double parse_nonnegative(std::string_view text) {
  std::size_t consumed = 0U;
  const double value = std::stod(std::string(text), &consumed);
  if (consumed != text.size() || !std::isfinite(value) || value < 0.0)
    throw std::runtime_error("invalid external performance metric");
  return value;
}

void append_external_metric(std::ostringstream& output,
                            std::string_view text,
                            std::string_view unit) {
  output << "{\"unit\":\"" << unit << "\",\"status\":";
  if (text == "null") {
    output << "\"unavailable\",\"value\":null}";
  } else {
    const double value = parse_nonnegative(text);
    output << "\"available\",\"value\":"
           << std::setprecision(std::numeric_limits<double>::max_digits10)
           << value << '}';
  }
}

int self_test();

int write_evidence(const std::filesystem::path& manifest,
                   const std::filesystem::path& destination) {
  std::ifstream stream(manifest);
  if (!stream)
    throw std::runtime_error("unable to read evidence manifest");
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << "{\"schema_version\":1,\"collection_boundary\":"
            "\"manual-release-evidence-only\",\"entries\":[";
  bool first = true;
  std::string line;
  while (std::getline(stream, line)) {
    if (line.empty())
      continue;
    const auto fields = split_tabs(line);
    if (fields.size() != 9U)
      throw std::runtime_error("evidence manifest field count differs");
    const auto case_bytes = read_file(fields[4U]);
    const auto artifact_bytes = read_file(fields[5U]);
    const auto input_bytes = read_file(fields[6U]);
    if (!hex64(fields[1U]) || !hex64(fields[2U]) ||
        !hex64(fields[3U]) || sha256(case_bytes) != fields[2U] ||
        sha256(artifact_bytes) != fields[3U])
      throw std::runtime_error("performance evidence SHA-256 differs");
    validate_input(input_bytes, artifact_bytes);
    if (!first)
      output << ',';
    first = false;
    output << "{\"case_id\":\"" << escape_json(fields[0U])
           << "\",\"command_sha256\":\"" << fields[1U]
           << "\",\"case_sha256\":\"" << fields[2U]
           << "\",\"artifact_sha256\":\"" << fields[3U]
           << "\",\"outer_wall_seconds\":";
    append_external_metric(output, fields[7U], "s");
    output << ",\"max_rss_kb\":";
    append_external_metric(output, fields[8U], "KiB");
    output << ",\"measurement\":" << input_bytes << '}';
  }
  if (first)
    throw std::runtime_error("evidence manifest is empty");
  output << "]}\n";
  const auto temporary =
      std::filesystem::path(destination.string() + ".tmp");
  try {
    std::ofstream destination_stream(
        temporary, std::ios::binary | std::ios::trunc);
    destination_stream << output.str();
    destination_stream.flush();
    if (!destination_stream)
      throw std::runtime_error("unable to stage evidence sidecar");
    destination_stream.close();
    std::filesystem::rename(temporary, destination);
  } catch (...) {
    std::error_code ignored;
    std::filesystem::remove(temporary, ignored);
    throw;
  }
  return 0;
}

std::string replace_once(std::string value, std::string_view before,
                         std::string_view after) {
  const auto offset = value.find(before);
  if (offset == std::string::npos)
    throw std::runtime_error("self-test mutation target is absent");
  value.replace(offset, before.size(), after);
  return value;
}

void expect_validation_failure(std::string_view input,
                               std::string_view artifact) {
  try {
    validate_input(input, artifact);
  } catch (const std::exception&) {
    return;
  }
  throw std::runtime_error("corrupt evidence passed validation");
}

void write_self_test_file(const std::filesystem::path& path,
                          std::string_view bytes) {
  std::ofstream stream(path, std::ios::binary | std::ios::trunc);
  stream.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  if (!stream)
    throw std::runtime_error("unable to write evidence self-test file");
}

int self_test() {
  if (sha256("") !=
          "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855" ||
      sha256("abc") !=
          "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
    throw std::runtime_error("SHA-256 self-test failed");

  hundun::diagnostics::PerformanceCorrectnessRecord correctness_record;
  correctness_record.passed = true;
  correctness_record.allocation_bytes_per_owned_cell = 1.0;
  correctness_record.peak_allocation_bytes_per_owned_cell = 2.0;
  correctness_record.repetitions = 2U;
  correctness_record.states = {
      {0U, "hundun-performance-state-fp-v1:0123456789abcdef"},
      {1U, "hundun-performance-state-fp-v1:0123456789abcdef"}};
  for (std::uint64_t repetition = 0U; repetition < 2U; ++repetition) {
    for (char phase : {'W', 'M'}) {
      const std::uint64_t steps = phase == 'W' ? 1U : 2U;
      for (std::uint64_t step = 0U; step < steps; ++step) {
        for (std::uint64_t slot = 0U; slot < 5U; ++slot) {
          correctness_record.work.push_back(
              {repetition, phase, step, slot, "converged", 1U, 2U, 3U,
               4U, UINT64_C(0x3ff0000000000000),
               UINT64_C(0x3fe0000000000000),
               UINT64_C(0x3fd0000000000000)});
        }
      }
    }
  }
  const auto correctness_summary =
      hundun::diagnostics::serialize_performance_correctness(
          correctness_record);
  const std::string artifact =
      "{\"schema_version\":1,\"correctness\":{\"passed\":true,"
      "\"summary\":\"" +
      escape_json(correctness_summary) +
      "\"},"
      "\"ranks\":1,\"process_grid\":[1,1,1],"
      "\"global_owned_cell_extents\":[8,8,4],"
      "\"per_rank_owned_cell_extents\":[8,8,4],"
      "\"measurement\":{\"warmup_steps\":1,\"measured_steps\":2,"
      "\"repetitions\":2},"
      "\"exact_counters\":{\"halo_payload_bytes\":{\"send\":2809856,"
      "\"receive\":2809856,\"pack\":2809856,\"unpack\":2809856},"
      "\"halo_messages\":{\"send\":3328,\"receive\":3328},"
      "\"logical_io_bytes\":{\"checkpoint\":0,\"diagnostics\":0}}}";
  const std::string runtime =
      "{\"completed_exchanges\":64,\"begin_calls\":64,\"wait_calls\":64,"
      "\"send_payload_bytes\":1404928,\"receive_payload_bytes\":1404928,"
      "\"pack_bytes\":1404928,\"unpack_bytes\":1404928,"
      "\"send_messages\":1664,\"receive_messages\":1664,"
      "\"successful_wait_seconds\":1.0}";
  const std::string pressure =
      "{\"completed_exchanges\":12,\"begin_calls\":12,\"wait_calls\":12,"
      "\"send_payload_bytes\":0,\"receive_payload_bytes\":0,"
      "\"pack_bytes\":0,\"unpack_bytes\":0,\"send_messages\":0,"
      "\"receive_messages\":0,\"successful_wait_seconds\":0.0}";
  const std::string combined =
      "{\"completed_exchanges\":76,\"begin_calls\":76,\"wait_calls\":76,"
      "\"send_payload_bytes\":1404928,\"receive_payload_bytes\":1404928,"
      "\"pack_bytes\":1404928,\"unpack_bytes\":1404928,"
      "\"send_messages\":1664,\"receive_messages\":1664,"
      "\"successful_wait_seconds\":1.0}";
  const std::string available =
      "{\"unit\":\"s\",\"status\":\"available\",\"value\":1.0}";
  const std::string bandwidth =
      "{\"unit\":\"byte/s\",\"status\":\"available\",\"value\":1404928.0}";
  const std::string unavailable_seconds =
      "{\"unit\":\"s\",\"status\":\"not_applicable\",\"value\":null}";
  const std::string unavailable_bandwidth =
      "{\"unit\":\"byte/s\",\"status\":\"not_applicable\",\"value\":null}";
  const std::string entry_zero =
      "{\"repetition\":0,"
      "\"relative_rank\":0,\"runtime\":" +
      runtime + ",\"pressure\":" + pressure + ",\"combined\":" + combined +
      "}";
  const std::string entry_one =
      replace_once(entry_zero, "\"repetition\":0", "\"repetition\":1");
  const std::string summaries_zero =
      "{\"repetition\":0,\"family\":\"runtime\","
      "\"critical_wait_seconds\":" +
      available + ",\"effective_bandwidth_bytes_per_second\":" + bandwidth +
      "},{\"repetition\":0,\"family\":\"pressure\","
      "\"critical_wait_seconds\":" +
      unavailable_seconds +
      ",\"effective_bandwidth_bytes_per_second\":" + unavailable_bandwidth +
      "},{\"repetition\":0,\"family\":\"combined\","
      "\"critical_wait_seconds\":" +
      available + ",\"effective_bandwidth_bytes_per_second\":" + bandwidth +
      "}";
  const std::string summaries_one =
      replace_once(replace_once(replace_once(
          summaries_zero, "\"repetition\":0", "\"repetition\":1"),
          "\"repetition\":0", "\"repetition\":1"),
          "\"repetition\":0", "\"repetition\":1");
  const std::string input =
      "{\"schema_version\":1,\"entries\":[" + entry_zero + "," +
      entry_one + "],\"halo_summaries\":[" + summaries_zero + "," +
      summaries_one +
      "],\"io\":{\"checkpoint\":{\"logical_bytes\":0,\"actual_bytes\":0,"
      "\"duration_seconds\":" +
      unavailable_seconds + ",\"throughput_bytes_per_second\":" +
      unavailable_bandwidth +
      "},\"diagnostics\":{\"logical_bytes\":0,\"actual_bytes\":0,"
      "\"duration_seconds\":" +
      unavailable_seconds + ",\"throughput_bytes_per_second\":" +
      unavailable_bandwidth + "}}}";
  validate_input(input, artifact);
  expect_validation_failure(
      input, replace_once(artifact, "\"exact_counters\"",
                           "\"wrong_counters\""));
  expect_validation_failure(
      input, replace_once(artifact, "\"ranks\":1", "\"ranks\":2"));
  expect_validation_failure(
      input, replace_once(artifact, "\"schema_version\":1",
                           "\"schema_version\":2"));
  expect_validation_failure(
      input, replace_once(artifact, "\"passed\":true",
                           "\"passed\":false"));
  expect_validation_failure(
      input, replace_once(
                 artifact,
                 "state=1:hundun-performance-state-fp-v1:0123456789abcdef",
                 "state=1:hundun-performance-state-fp-v1:1123456789abcdef"));
  expect_validation_failure(
      input, replace_once(artifact, "work=1:W:0:0:converged:1:",
                           "work=1:W:0:0:converged:9:"));
  expect_validation_failure(
      replace_once(input, "\"completed_exchanges\":64",
                   "\"completed_exchanges\":65"),
      artifact);
  expect_validation_failure(
      replace_once(input, "\"relative_rank\":0",
                   "\"relative_rank\":18446744073709551615"),
      artifact);
  expect_validation_failure(
      replace_once(input, "\"actual_bytes\":0",
                   "\"actual_bytes\":1"),
      artifact);
  bool overflow_rejected = false;
  try {
    static_cast<void>(
        checked_add(std::numeric_limits<std::uint64_t>::max(), 1U));
  } catch (const std::runtime_error&) {
    overflow_rejected = true;
  }
  if (!overflow_rejected)
    throw std::runtime_error("overflow self-test unexpectedly succeeded");

  std::filesystem::path directory;
  for (unsigned attempt = 0U; attempt < 100U; ++attempt) {
    directory = std::filesystem::temp_directory_path() /
                ("hundun-task25-evidence-self-" +
                 std::to_string(attempt));
    std::error_code create_error;
    if (std::filesystem::create_directory(directory, create_error))
      break;
    directory.clear();
  }
  if (directory.empty())
    throw std::runtime_error("unable to create evidence self-test directory");
  try {
    const auto case_path = directory / "case.json";
    const auto artifact_path = directory / "artifact.json";
    const auto input_path = directory / "input.json";
    const auto manifest_path = directory / "manifest.tsv";
    const auto destination = directory / "evidence.json";
    write_self_test_file(case_path, "{}");
    write_self_test_file(artifact_path, artifact);
    write_self_test_file(input_path, input);
    write_self_test_file(destination, "preserved");
    const std::string wrong_sha(64U, '0');
    const std::string command_sha(64U, '1');
    const std::string manifest =
        "self-test\t" + command_sha + "\t" + sha256("{}") + "\t" +
        wrong_sha + "\t" + case_path.string() + "\t" +
        artifact_path.string() + "\t" + input_path.string() +
        "\tnull\tnull\n";
    write_self_test_file(manifest_path, manifest);
    bool sha_rejected = false;
    try {
      static_cast<void>(write_evidence(manifest_path, destination));
    } catch (const std::runtime_error&) {
      sha_rejected = true;
    }
    if (!sha_rejected)
      throw std::runtime_error("SHA mismatch self-test unexpectedly succeeded");
    if (read_file(destination) != "preserved" ||
        std::filesystem::exists(
            std::filesystem::path(destination.string() + ".tmp")))
      throw std::runtime_error("failed evidence write changed destination");
  } catch (...) {
    std::error_code cleanup_error;
    std::filesystem::remove_all(directory, cleanup_error);
    throw;
  }
  std::error_code cleanup_error;
  std::filesystem::remove_all(directory, cleanup_error);
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string_view(argv[1]) == "--self-test")
      return self_test();
    if (argc != 3) {
      std::cerr << "usage: task25_performance_evidence "
                   "<manifest.tsv> <performance-evidence.v1.json>\n";
      return 2;
    }
    return write_evidence(argv[1], argv[2]);
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
