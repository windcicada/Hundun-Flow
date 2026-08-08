// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_performance_correctness.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>

namespace hundun::diagnostics {
namespace {

constexpr std::string_view kPrefix = "hundun-performance-correctness-v1";
constexpr std::string_view kStatePrefix =
    "hundun-performance-state-fp-v1:";

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

double from_bits(std::uint64_t value) noexcept {
  double result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

std::string hex16(std::uint64_t value) {
  constexpr char digits[] = "0123456789abcdef";
  std::string result(16U, '0');
  for (std::size_t index = 0; index < result.size(); ++index) {
    const auto shift = static_cast<unsigned>((15U - index) * 4U);
    result[index] = digits[(value >> shift) & UINT64_C(0xf)];
  }
  return result;
}

std::uint64_t parse_u64(std::string_view text) {
  if (text.empty() || (text.size() > 1U && text.front() == '0')) {
    throw std::invalid_argument("performance integer is not minimal base-10");
  }
  std::uint64_t result = 0U;
  for (char character : text) {
    if (character < '0' || character > '9') {
      throw std::invalid_argument("performance integer is not base-10");
    }
    const auto digit = static_cast<std::uint64_t>(character - '0');
    if (result > (std::numeric_limits<std::uint64_t>::max() - digit) / 10U) {
      throw std::invalid_argument("performance integer overflows");
    }
    result = result * 10U + digit;
  }
  return result;
}

std::uint64_t parse_hex16(std::string_view text) {
  if (text.size() != 16U) {
    throw std::invalid_argument("performance FP64 bits require 16 hex digits");
  }
  std::uint64_t result = 0U;
  for (char character : text) {
    std::uint64_t digit = 0U;
    if (character >= '0' && character <= '9') {
      digit = static_cast<std::uint64_t>(character - '0');
    } else if (character >= 'a' && character <= 'f') {
      digit = static_cast<std::uint64_t>(character - 'a' + 10);
    } else {
      throw std::invalid_argument(
          "performance FP64 bits require lowercase hexadecimal");
    }
    result = result * 16U + digit;
  }
  return result;
}

bool valid_termination(std::string_view value) noexcept {
  constexpr std::array<std::string_view, 7> valid{
      "converged",          "zero_right_hand_side", "maximum_iterations",
      "numerical_breakdown", "non_finite_value",     "invalid_control",
      "collective_failure"};
  return std::find(valid.begin(), valid.end(), value) != valid.end();
}

bool successful_termination(std::string_view value) noexcept {
  return value == "converged" || value == "zero_right_hand_side";
}

bool valid_residual_bits(std::uint64_t value) noexcept {
  const double residual = from_bits(value);
  return std::isfinite(residual) && residual >= 0.0;
}

bool same_work_evidence(const PerformanceWorkRecord& left,
                        const PerformanceWorkRecord& right) noexcept {
  return left.phase == right.phase &&
         left.relative_step == right.relative_step &&
         left.slot == right.slot &&
         left.termination == right.termination &&
         left.iterations == right.iterations &&
         left.matvec == right.matvec &&
         left.preconditioner == right.preconditioner &&
         left.reduction == right.reduction &&
         left.initial_residual_bits == right.initial_residual_bits &&
         left.recursive_residual_bits == right.recursive_residual_bits &&
         left.independent_final_residual_bits ==
             right.independent_final_residual_bits;
}

std::vector<std::string_view> split(std::string_view input, char separator) {
  std::vector<std::string_view> result;
  std::size_t begin = 0U;
  for (;;) {
    const auto end = input.find(separator, begin);
    result.push_back(input.substr(
        begin, end == std::string_view::npos ? input.size() - begin
                                             : end - begin));
    if (end == std::string_view::npos) {
      return result;
    }
    begin = end + 1U;
  }
}

void validate_state_fingerprint(std::string_view value) {
  if (value.size() != kStatePrefix.size() + 16U ||
      value.substr(0U, kStatePrefix.size()) != kStatePrefix) {
    throw std::invalid_argument("invalid performance state fingerprint");
  }
  static_cast<void>(parse_hex16(value.substr(kStatePrefix.size())));
}

void validate(const PerformanceCorrectnessRecord& record) {
  if (!std::isfinite(record.allocation_bytes_per_owned_cell) ||
      record.allocation_bytes_per_owned_cell < 0.0 ||
      !std::isfinite(record.peak_allocation_bytes_per_owned_cell) ||
      record.peak_allocation_bytes_per_owned_cell < 0.0 ||
      record.peak_allocation_bytes_per_owned_cell <
          record.allocation_bytes_per_owned_cell ||
      record.repetitions == 0U ||
      record.states.size() !=
          static_cast<std::size_t>(record.repetitions)) {
    throw std::invalid_argument("invalid performance correctness header");
  }
  for (std::size_t index = 0U; index < record.states.size(); ++index) {
    if (record.states[index].first != static_cast<std::uint64_t>(index)) {
      throw std::invalid_argument(
          "performance states must cover repetitions in order");
    }
    validate_state_fingerprint(record.states[index].second);
    if (record.passed && index != 0U &&
        record.states[index].second != record.states.front().second) {
      throw std::invalid_argument(
          "successful performance states differ between repetitions");
    }
  }

  std::optional<std::array<std::uint64_t, 2>> reference_steps;
  std::size_t index = 0U;
  for (std::uint64_t repetition = 0U; repetition < record.repetitions;
       ++repetition) {
    std::array<std::uint64_t, 2> steps{};
    for (int phase_index = 0; phase_index < 2; ++phase_index) {
      const char phase = phase_index == 0 ? 'W' : 'M';
      std::uint64_t relative_step = 0U;
      while (index < record.work.size() &&
             record.work[index].repetition == repetition &&
             record.work[index].phase == phase) {
        if (record.work.size() - index < 5U) {
          throw std::invalid_argument(
              "performance work must contain complete five-report groups");
        }
        for (std::uint64_t slot = 0U; slot < 5U; ++slot) {
          const auto& item = record.work[index++];
          if (item.repetition != repetition || item.phase != phase ||
              item.relative_step != relative_step || item.slot != slot ||
              !valid_termination(item.termination)) {
            throw std::invalid_argument(
                "performance work group key or slot is invalid");
          }
        }
        ++relative_step;
      }
      steps[static_cast<std::size_t>(phase_index)] = relative_step;
    }
    if (index < record.work.size() &&
        record.work[index].repetition == repetition) {
      throw std::invalid_argument(
          "performance work phase or step ordering is invalid");
    }
    if (!reference_steps) {
      reference_steps = steps;
    } else if (*reference_steps != steps) {
      throw std::invalid_argument(
          "performance repetitions must have identical step coverage");
    }
  }
  if (index != record.work.size()) {
    throw std::invalid_argument(
        "performance work repetition ordering is invalid");
  }
  for (const auto& item : record.work) {
    if (item.repetition >= record.repetitions ||
        (item.phase != 'W' && item.phase != 'M') || item.slot >= 5U ||
        !valid_termination(item.termination) ||
        (record.passed && !successful_termination(item.termination)) ||
        !valid_residual_bits(item.initial_residual_bits) ||
        !valid_residual_bits(item.recursive_residual_bits) ||
        !valid_residual_bits(item.independent_final_residual_bits)) {
      throw std::invalid_argument("invalid performance work item");
    }
  }
  if (record.work.empty() || record.work.size() % 5U != 0U) {
    throw std::invalid_argument(
        "performance work must contain complete five-report groups");
  }
  if (record.passed) {
    const auto per_repetition =
        record.work.size() /
        static_cast<std::size_t>(record.repetitions);
    for (std::uint64_t repetition = 1U;
         repetition < record.repetitions; ++repetition) {
      const auto offset =
          static_cast<std::size_t>(repetition) * per_repetition;
      for (std::size_t work_index = 0U;
           work_index < per_repetition; ++work_index) {
        if (!same_work_evidence(record.work[work_index],
                                record.work[offset + work_index])) {
          throw std::invalid_argument(
              "successful performance work differs between repetitions");
        }
      }
    }
  }
}

}  // namespace

std::string serialize_performance_correctness(
    const PerformanceCorrectnessRecord& record) {
  validate(record);
  std::string result(kPrefix);
  result += ";passed=";
  result += record.passed ? '1' : '0';
  result += ";allocation-bytes-per-owned-cell=" +
            hex16(bits(record.allocation_bytes_per_owned_cell));
  result += ";peak-allocation-bytes-per-owned-cell=" +
            hex16(bits(record.peak_allocation_bytes_per_owned_cell));
  result += ";repetitions=" + std::to_string(record.repetitions);
  for (const auto& state : record.states) {
    result += ";state=" + std::to_string(state.first) + ':' + state.second;
  }
  for (const auto& item : record.work) {
    result += ";work=" + std::to_string(item.repetition);
    result.push_back(':');
    result.push_back(item.phase);
    result += ':' + std::to_string(item.relative_step) + ':' +
              std::to_string(item.slot) + ':' + item.termination + ':' +
              std::to_string(item.iterations) + ':' +
              std::to_string(item.matvec) + ':' +
              std::to_string(item.preconditioner) + ':' +
              std::to_string(item.reduction) + ':' +
              hex16(item.initial_residual_bits) + ':' +
              hex16(item.recursive_residual_bits) + ':' +
              hex16(item.independent_final_residual_bits);
  }
  return result;
}

PerformanceCorrectnessRecord parse_performance_correctness(
    std::string_view encoded) {
  if (encoded.empty() || encoded.size() > UINT64_C(64) * 1024U * 1024U) {
    throw std::invalid_argument("performance correctness record size invalid");
  }
  const auto fields = split(encoded, ';');
  if (fields.size() < 7U || fields.front() != kPrefix) {
    throw std::invalid_argument("invalid performance correctness prefix");
  }
  PerformanceCorrectnessRecord result;
  std::size_t index = 1U;
  const auto take = [&](std::string_view key) {
    if (index >= fields.size() ||
        fields[index].substr(0U, key.size()) != key) {
      throw std::invalid_argument("performance correctness field missing");
    }
    return fields[index++].substr(key.size());
  };
  const auto passed = take("passed=");
  if (passed != "0" && passed != "1") {
    throw std::invalid_argument("invalid performance correctness passed flag");
  }
  result.passed = passed == "1";
  result.allocation_bytes_per_owned_cell = from_bits(
      parse_hex16(take("allocation-bytes-per-owned-cell=")));
  result.peak_allocation_bytes_per_owned_cell =
      from_bits(parse_hex16(take("peak-allocation-bytes-per-owned-cell=")));
  result.repetitions = parse_u64(take("repetitions="));

  while (index < fields.size() &&
         fields[index].substr(0U, 6U) == "state=") {
    const auto value = fields[index++].substr(6U);
    const auto separator = value.find(':');
    if (separator == std::string_view::npos) {
      throw std::invalid_argument("invalid performance state item");
    }
    result.states.emplace_back(parse_u64(value.substr(0U, separator)),
                               std::string(value.substr(separator + 1U)));
  }
  while (index < fields.size()) {
    if (fields[index].substr(0U, 5U) != "work=") {
      throw std::invalid_argument("unknown performance correctness field");
    }
    const auto parts = split(fields[index++].substr(5U), ':');
    if (parts.size() != 12U || parts[1].size() != 1U) {
      throw std::invalid_argument("invalid performance work item");
    }
    PerformanceWorkRecord item;
    item.repetition = parse_u64(parts[0]);
    item.phase = parts[1].front();
    item.relative_step = parse_u64(parts[2]);
    item.slot = parse_u64(parts[3]);
    item.termination = std::string(parts[4]);
    item.iterations = parse_u64(parts[5]);
    item.matvec = parse_u64(parts[6]);
    item.preconditioner = parse_u64(parts[7]);
    item.reduction = parse_u64(parts[8]);
    item.initial_residual_bits = parse_hex16(parts[9]);
    item.recursive_residual_bits = parse_hex16(parts[10]);
    item.independent_final_residual_bits = parse_hex16(parts[11]);
    result.work.push_back(std::move(item));
  }
  validate(result);
  return result;
}

void validate_performance_correctness_coverage(
    const PerformanceCorrectnessRecord& record, std::uint64_t warmup_steps,
    std::uint64_t measured_steps, std::uint64_t repetitions) {
  validate(record);
  if (record.repetitions != repetitions) {
    throw std::invalid_argument(
        "performance correctness repetition metadata mismatch");
  }
  if (warmup_steps >
          std::numeric_limits<std::uint64_t>::max() - measured_steps ||
      warmup_steps + measured_steps >
          std::numeric_limits<std::uint64_t>::max() / 5U ||
      (repetitions != 0U &&
       5U * (warmup_steps + measured_steps) >
           std::numeric_limits<std::uint64_t>::max() / repetitions)) {
    throw std::invalid_argument(
        "performance correctness coverage overflows");
  }
  const auto expected =
      repetitions * 5U * (warmup_steps + measured_steps);
  if (expected != static_cast<std::uint64_t>(record.work.size())) {
    throw std::invalid_argument(
        "performance correctness work coverage mismatch");
  }
  std::size_t index = 0U;
  for (std::uint64_t repetition = 0U; repetition < repetitions;
       ++repetition) {
    for (const auto& phase :
         std::array<std::pair<char, std::uint64_t>, 2>{
             std::pair{'W', warmup_steps},
             std::pair{'M', measured_steps}}) {
      for (std::uint64_t step = 0U; step < phase.second; ++step) {
        for (std::uint64_t slot = 0U; slot < 5U; ++slot) {
          const auto& item = record.work[index++];
          if (item.repetition != repetition ||
              item.phase != phase.first || item.relative_step != step ||
              item.slot != slot) {
            throw std::invalid_argument(
                "performance correctness work metadata mismatch");
          }
        }
      }
    }
  }
}

}  // namespace hundun::diagnostics
