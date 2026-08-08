// SPDX-License-Identifier: Apache-2.0

#include "app_cli_options_detail.hpp"

#include "hundun/rt_error.hpp"

#include <string_view>

namespace hundun::application {
namespace {

[[noreturn]] void throw_usage() { throw runtime::Error(kCliUsage); }

bool is_case_path(std::string_view argument) {
  return !argument.empty() && argument.front() != '-';
}

}  // namespace

CliOptions parse_cli(int argc, char** argv) {
  if (argv == nullptr || argc < 2) {
    throw_usage();
  }

  if (argc == 2 && argv[1] != nullptr) {
    const std::string_view argument{argv[1]};
    if (argument == "--version") {
      CliOptions options;
      options.show_version = true;
      return options;
    }
    if (is_case_path(argument)) {
      CliOptions options;
      options.case_path = argument;
      return options;
    }
    throw_usage();
  }

  if (argc == 3 && argv[1] != nullptr && argv[2] != nullptr) {
    const std::string_view case_path{argv[1]};
    const std::string_view mode{argv[2]};
    if (!is_case_path(case_path)) {
      throw_usage();
    }

    CliOptions options;
    options.case_path = case_path;
    if (mode == "--validate") {
      options.validate_only = true;
      return options;
    }
    if (mode == "--print-resolved") {
      options.print_resolved = true;
      return options;
    }
  }

  throw_usage();
}

}  // namespace hundun::application
