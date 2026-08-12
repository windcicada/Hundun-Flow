// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/cfg_resolved_case_v4.hpp"

#include <filesystem>
#include <string>

namespace hundun::config {

ResolvedReactingCaseV4
load_resolved_reacting_case_v4(const std::filesystem::path &path);

std::string
to_resolved_reacting_json_v4(const ResolvedReactingCaseV4 &resolved_case);

} // namespace hundun::config
