// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/cfg_case_config.hpp"

#include <filesystem>
#include <string>

namespace hundun::config {

CaseConfig load_case_config(const std::filesystem::path& path);
void validate_case_config(const CaseConfig& config);
std::string to_resolved_json(const CaseConfig& config);

}  // namespace hundun::config
