// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#error "structured diagnostics test access is test-build only"
#endif

#include <cstdint>
#include <string_view>

namespace hundun::diagnostics::test {

std::uint64_t crc64_ecma(std::string_view bytes) noexcept;

} // namespace hundun::diagnostics::test
