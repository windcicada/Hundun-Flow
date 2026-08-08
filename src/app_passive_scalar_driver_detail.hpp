// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "app_cli_options_detail.hpp"
#include "hundun/cfg_case_config.hpp"
#include "hundun/rt_mpi_context.hpp"

#include <filesystem>
#include <string_view>

namespace hundun::application {

std::filesystem::path establish_authoritative_case_root(
    const runtime::MpiContext& context,
    const std::filesystem::path& local_case_path);

int run_passive_scalar_case(const CliOptions& options, runtime::MpiContext& context,
                    const config::CaseConfig& config,
                    const std::filesystem::path& authoritative_case_root);

}  // namespace hundun::application
