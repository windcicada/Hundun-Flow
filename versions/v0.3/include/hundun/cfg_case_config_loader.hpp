// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/cfg_case_config.hpp"

#include <filesystem>
#include <string>

namespace hundun::config {

CaseConfig load_case_config(const std::filesystem::path& path);
void validate_case_config(const CaseConfig& config);
std::string to_resolved_json(const CaseConfig& config);

}  // namespace hundun::config
