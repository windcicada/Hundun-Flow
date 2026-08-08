// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/cfg_resolved_case_v3.hpp"

#include <filesystem>
#include <string>

namespace hundun::config::detail {

ResolvedCaseV3
parse_resolved_case_v3_json(std::string contents,
                            const std::filesystem::path &authoritative_path);

} // namespace hundun::config::detail
