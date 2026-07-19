// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/config/case_config.hpp"

#include <filesystem>
#include <string_view>

namespace hundun::config::detail {

CaseConfig load_case_config_from_json(
    std::string_view contents, const std::filesystem::path& case_path);

}  // namespace hundun::config::detail
