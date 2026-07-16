// SPDX-License-Identifier: Apache-2.0

#include "applications/hundun/cli_options.hpp"
#include "runtime/include/hundun/runtime/error.hpp"
#include "tests/support/test_main.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace {

using hundun::application::CliOptions;
using hundun::application::kCliUsage;
using hundun::application::parse_cli;

constexpr std::string_view kExpectedUsage =
    "usage: hundun <case.json> [--validate|--print-resolved] | hundun --version";

CliOptions parse(std::vector<std::string> arguments) {
  std::vector<char*> argv;
  argv.reserve(arguments.size());
  for (auto& argument : arguments) {
    argv.push_back(argument.data());
  }
  return parse_cli(static_cast<int>(argv.size()), argv.data());
}

void check_defaults(const CliOptions& options) {
  HUNDUN_CHECK(!options.validate_only);
  HUNDUN_CHECK(!options.print_resolved);
  HUNDUN_CHECK(!options.show_version);
}

void test_accepted_forms() {
  const CliOptions run = parse({"hundun", "case.json"});
  HUNDUN_CHECK(run.case_path == "case.json");
  check_defaults(run);

  const CliOptions validate =
      parse({"hundun", "nested/case.json", "--validate"});
  HUNDUN_CHECK(validate.case_path == "nested/case.json");
  HUNDUN_CHECK(validate.validate_only);
  HUNDUN_CHECK(!validate.print_resolved);
  HUNDUN_CHECK(!validate.show_version);

  const CliOptions resolved =
      parse({"hundun", "case.json", "--print-resolved"});
  HUNDUN_CHECK(resolved.case_path == "case.json");
  HUNDUN_CHECK(!resolved.validate_only);
  HUNDUN_CHECK(resolved.print_resolved);
  HUNDUN_CHECK(!resolved.show_version);

  const CliOptions version = parse({"hundun", "--version"});
  HUNDUN_CHECK(version.case_path.empty());
  HUNDUN_CHECK(!version.validate_only);
  HUNDUN_CHECK(!version.print_resolved);
  HUNDUN_CHECK(version.show_version);
}

void expect_usage_error(std::vector<std::string> arguments) {
  bool threw = false;
  try {
    (void)parse(std::move(arguments));
  } catch (const hundun::runtime::Error& error) {
    threw = true;
    HUNDUN_CHECK(std::string_view(error.what()) == kExpectedUsage);
  }
  HUNDUN_CHECK(threw);
}

void test_usage_contract() {
  HUNDUN_CHECK(std::string_view(kCliUsage) == kExpectedUsage);
}

void test_rejected_forms() {
  expect_usage_error({"hundun"});
  expect_usage_error({"hundun", ""});
  expect_usage_error({"hundun", "--unknown"});
  expect_usage_error({"hundun", "case.json", "--unknown"});
  expect_usage_error({"hundun", "first.json", "second.json"});
  expect_usage_error({"hundun", "case.json", "--validate",
                      "--print-resolved"});
  expect_usage_error({"hundun", "case.json", "--print-resolved",
                      "--validate"});
  expect_usage_error({"hundun", "--validate", "case.json"});
  expect_usage_error({"hundun", "--print-resolved", "case.json"});
  expect_usage_error({"hundun", "--validate"});
  expect_usage_error({"hundun", "--print-resolved"});
  expect_usage_error({"hundun", "case.json", "--validate", "--validate"});
  expect_usage_error(
      {"hundun", "case.json", "--print-resolved", "--print-resolved"});
  expect_usage_error({"hundun", "--version", "case.json"});
  expect_usage_error({"hundun", "case.json", "--version"});
  expect_usage_error({"hundun", "--version", "--validate"});
  expect_usage_error({"hundun", "--version", "--version"});
  expect_usage_error({"hundun", "case.json", ""});
  expect_usage_error({"hundun", "case.json", "--validate", "extra.json"});
}

}  // namespace

int main() {
  return hundun::test::run([] {
    test_usage_contract();
    test_accepted_forms();
    test_rejected_forms();
  });
}
