// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "applications/hundun/cli_options.hpp"

#include "hundun/config/case_config.hpp"
#include "hundun/runtime/mpi_context.hpp"

#include <filesystem>
#include <string_view>

namespace hundun::application {

inline constexpr std::string_view kHundunVersion =
    "HUNDUN-FLOW 0.0.0-stage1";

std::filesystem::path establish_authoritative_case_root(
    const runtime::MpiContext& context,
    const std::filesystem::path& local_case_path);

int run_stage1_case(const CliOptions& options, runtime::MpiContext& context,
                    const config::CaseConfig& config,
                    const std::filesystem::path& authoritative_case_root);

}  // namespace hundun::application
