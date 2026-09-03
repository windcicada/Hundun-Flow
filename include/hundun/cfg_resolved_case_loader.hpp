// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/cfg_case_config_loader.hpp"
#include "hundun/cfg_resolved_case.hpp"

#include <filesystem>
#include <string>

namespace hundun::config {

ResolvedCase load_resolved_case(const std::filesystem::path& path);
std::string to_resolved_json(const ResolvedCase& resolved_case);

}  // namespace hundun::config
