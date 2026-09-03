// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "app_cli_options_detail.hpp"

#include "hundun/cfg_resolved_case.hpp"
#include "hundun/rt_mpi_context.hpp"

#include <filesystem>

namespace hundun::application {

int run_flow_case(
    const CliOptions& options, runtime::MpiContext& context,
    const config::FlowCaseConfig& config,
    const std::filesystem::path& authoritative_case_root);

}  // namespace hundun::application
