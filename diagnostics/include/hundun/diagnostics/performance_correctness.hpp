// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace hundun::diagnostics {

struct PerformanceWorkRecord final {
  std::uint64_t repetition{};
  char phase{'W'};
  std::uint64_t relative_step{};
  std::uint64_t slot{};
  std::string termination;
  std::uint64_t iterations{};
  std::uint64_t matvec{};
  std::uint64_t preconditioner{};
  std::uint64_t reduction{};
  std::uint64_t initial_residual_bits{};
  std::uint64_t recursive_residual_bits{};
  std::uint64_t independent_final_residual_bits{};
};

struct PerformanceCorrectnessRecord final {
  bool passed{};
  double allocation_bytes_per_owned_cell{};
  double peak_allocation_bytes_per_owned_cell{};
  std::uint64_t repetitions{};
  std::vector<std::pair<std::uint64_t, std::string>> states;
  std::vector<PerformanceWorkRecord> work;
};

std::string serialize_performance_correctness(
    const PerformanceCorrectnessRecord& record);

PerformanceCorrectnessRecord parse_performance_correctness(
    std::string_view encoded);

void validate_performance_correctness_coverage(
    const PerformanceCorrectnessRecord& record, std::uint64_t warmup_steps,
    std::uint64_t measured_steps, std::uint64_t repetitions);

}  // namespace hundun::diagnostics
