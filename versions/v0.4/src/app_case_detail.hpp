// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_case.hpp"

#include <cstddef>
#include <cstdint>

namespace hundun::v04::detail {

inline constexpr std::size_t kMaxJsonBytes = 1024U * 1024U;
inline constexpr std::size_t kMaxReferencedFiles = 256U;
inline constexpr std::size_t kMaxRelativePathBytes = 255U;
inline constexpr std::uint64_t kMaxReferencedFileBytes = 64U * 1024U * 1024U;
inline constexpr std::size_t kMaxWireBytes = 128U * 1024U;

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
using FileOpenObserver = void (*)(int rank);
void set_file_open_observer_for_test(FileOpenObserver observer) noexcept;
int last_lowest_failing_rank_for_test() noexcept;
#endif

}  // namespace hundun::v04::detail
