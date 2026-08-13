// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_mesh.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::v04::detail {

inline constexpr std::uint64_t kMaxStlBytes = UINT64_C(1073741824);
inline constexpr std::uint64_t kMaxStlTriangles = UINT64_C(20000000);

enum StlScanDetail : std::uint32_t {
  stl_detail_none = 0U,
  stl_detail_geometry = 1U,
  stl_detail_patch = 2U,
  stl_detail_triangle = 3U,
  stl_detail_open_surface = 4U,
  stl_detail_path = 5U,
  stl_detail_open = 6U,
  stl_detail_read = 7U,
  stl_detail_size = 8U,
  stl_detail_syntax = 9U,
  stl_detail_collective = 10U,
  stl_detail_region = 11U,
  stl_detail_allocation = 12U,
  stl_detail_budget = 13U,
  stl_detail_bin_references = 14U,
};

Status parse_stl_bytes(Span<const std::uint8_t> bytes,
                       std::vector<TriangleInput>& triangles) noexcept;

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
void reset_stl_open_count_for_test() noexcept;
std::uint64_t stl_open_count_for_test() noexcept;
void fail_thread_launch_after_for_test(std::size_t successful_launches) noexcept;
void reset_thread_launch_failure_for_test() noexcept;
#endif

struct LineChunk {
  std::size_t begin{};
  std::size_t end{};
};

LineChunk fixed_line_chunk(std::size_t line_count, std::size_t worker,
                           std::size_t worker_count) noexcept;

}  // namespace hundun::v04::detail
