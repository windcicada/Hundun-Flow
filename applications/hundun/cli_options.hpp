// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <filesystem>

namespace hundun::application {

inline constexpr char kCliUsage[] =
    "usage: hundun <case.json> [--validate|--print-resolved] | hundun "
    "--version";

struct CliOptions {
  std::filesystem::path case_path;
  bool validate_only{false};
  bool print_resolved{false};
  bool show_version{false};
};

CliOptions parse_cli(int argc, char** argv);

}  // namespace hundun::application
