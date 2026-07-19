// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/config/case_config_loader.hpp"
#include "hundun/config/resolved_case.hpp"

#include <filesystem>
#include <string>

namespace hundun::config {

ResolvedCase load_resolved_case(const std::filesystem::path& path);
std::string to_resolved_json(const ResolvedCase& resolved_case);

}  // namespace hundun::config
