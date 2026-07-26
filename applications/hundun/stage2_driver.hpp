// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "applications/hundun/cli_options.hpp"

#include "hundun/config/resolved_case.hpp"
#include "hundun/runtime/mpi_context.hpp"

#include <filesystem>

namespace hundun::application {

int run_stage2_case(
    const CliOptions& options, runtime::MpiContext& context,
    const config::FlowCaseConfig& config,
    const std::filesystem::path& authoritative_case_root);

}  // namespace hundun::application
