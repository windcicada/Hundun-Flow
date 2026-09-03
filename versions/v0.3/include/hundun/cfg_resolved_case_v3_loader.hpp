// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/cfg_resolved_case_v3.hpp"

#include <filesystem>
#include <string>

namespace hundun::config {

ResolvedCaseV3 load_resolved_case_v3(const std::filesystem::path &path);

std::string to_resolved_json_v3(const ResolvedCaseV3 &resolved_case);

} // namespace hundun::config
