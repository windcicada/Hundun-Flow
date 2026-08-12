// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/cfg_resolved_case_v4.hpp"

#include <filesystem>
#include <string>

namespace hundun::config::detail {

ResolvedReactingCaseV4 parse_resolved_reacting_case_v4_json(
    std::string contents, const std::filesystem::path &authoritative_path);

} // namespace hundun::config::detail
